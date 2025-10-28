#!/usr/bin/env python3
"""
Collaborative Analysis Server for RNASEQ-MINI

This server integrates the collaboration manager, notification system,
and Jupyter integration to provide a complete real-time collaborative
analysis environment.
"""

import asyncio
import argparse
import threading
import time
from pathlib import Path
from scripts.collaboration_manager import CollaborationManager
from scripts.notification_manager import NotificationManager
from scripts.jupyter_integration import JupyterIntegration


class CollaborativeAnalysisServer:
    """Main server for collaborative RNA-seq analysis."""

    def __init__(self, host: str = "localhost", port: int = 8765,
                 enable_notifications: bool = True, enable_jupyter: bool = True):
        self.host = host
        self.port = port
        self.enable_notifications = enable_notifications
        self.enable_jupyter = enable_jupyter

        # Initialize components
        self.collaboration_manager = CollaborationManager()
        self.jupyter_integration = JupyterIntegration(self.collaboration_manager) if enable_jupyter else None

        if enable_notifications:
            self.notification_manager = NotificationManager()
            self.collaboration_manager.set_notification_manager(self.notification_manager)
        else:
            self.notification_manager = None

        # Track running tasks
        self.server_task = None
        self.monitoring_task = None

    async def start_server(self):
        """Start the collaborative server."""
        print("🚀 Starting RNASEQ-MINI Collaborative Analysis Server")
        print(f"🛜 WebSocket server: ws://{self.host}:{self.port}")
        print(f"📢 Notifications: {'Enabled' if self.enable_notifications else 'Disabled'}")
        print(f"📓 Jupyter integration: {'Enabled' if self.enable_jupyter else 'Disabled'}")

        # Start collaboration WebSocket server
        server = await self.collaboration_manager.start_collaboration_server(self.host, self.port)

        # Start monitoring task
        self.monitoring_task = asyncio.create_task(self._monitoring_loop())

        try:
            await server.wait_closed()
        except KeyboardInterrupt:
            print("\n🛑 Shutting down server...")
        finally:
            # Cancel monitoring task
            if self.monitoring_task:
                self.monitoring_task.cancel()

    async def _monitoring_loop(self):
        """Background monitoring loop for server health."""
        while True:
            try:
                await asyncio.sleep(60)  # Check every minute

                # Log server status
                active_sessions = len(self.collaboration_manager.active_sessions)
                total_connections = sum(len(conns) for conns in self.collaboration_manager.websocket_connections.values())

                print(f"📊 Server Status - Sessions: {active_sessions}, Connections: {total_connections}")

                # Check for inactive sessions to clean up
                await self._cleanup_inactive_sessions()

            except asyncio.CancelledError:
                break
            except Exception as e:
                print(f"❌ Monitoring error: {e}")

    async def _cleanup_inactive_sessions(self):
        """Clean up inactive sessions."""
        # This would implement session cleanup logic
        # For now, just log that we're checking
        pass

    def create_session(self, project_name: str, config_file: str = "config/params.yaml") -> str:
        """Create a new collaborative session."""
        return self.collaboration_manager.create_session(
            name=f"Collaborative: {project_name}",
            description=f"Shared RNA-seq analysis for {project_name}",
            created_by="server",
            project_config=self._load_config(config_file)
        )

    def _load_config(self, config_file: str) -> dict:
        """Load project configuration."""
        try:
            import yaml
            with open(config_file) as f:
                return yaml.safe_load(f)
        except Exception as e:
            print(f"❌ Failed to load config {config_file}: {e}")
            return {}

    def setup_notifications(self):
        """Interactive setup for notification channels."""
        if not self.notification_manager:
            print("❌ Notifications not enabled")
            return

        from scripts.notification_manager import setup_notifications
        setup_notifications()

    def create_collaborative_notebook(self, session_id: str, template: str = "standard") -> str:
        """Create a Jupyter notebook for collaborative analysis."""
        if not self.jupyter_integration:
            print("❌ Jupyter integration not enabled")
            return None

        return self.jupyter_integration.create_analysis_notebook(session_id, template)

    def list_active_sessions(self):
        """List all active collaborative sessions."""
        sessions = []
        for session_id, session in self.collaboration_manager.active_sessions.items():
            sessions.append({
                "id": session_id,
                "name": session.name,
                "participants": len(session.participants),
                "current_step": session.current_step,
                "created": session.created_at.isoformat()
            })

        return sessions

    def get_session_info(self, session_id: str) -> dict:
        """Get detailed information about a session."""
        session = self.collaboration_manager.db.get_session(session_id)
        if not session:
            return None

        # Get recent events
        events = self.collaboration_manager.db.get_session_events(session_id, 10)

        return {
            "session": {
                "id": session.id,
                "name": session.name,
                "description": session.description,
                "created_by": session.created_by,
                "created_at": session.created_at.isoformat(),
                "status": session.status,
                "current_step": session.current_step,
                "participants": {uid: {
                    "name": user.name,
                    "role": user.role,
                    "joined_at": user.joined_at.isoformat(),
                    "last_active": user.last_active.isoformat()
                } for uid, user in session.participants.items()}
            },
            "recent_events": [{
                "id": event.id,
                "type": event.event_type,
                "description": event.description,
                "user": next((p.name for p in session.participants.values() if p.id == event.user_id), "Unknown"),
                "timestamp": event.timestamp.isoformat(),
                "step": event.step
            } for event in events]
        }


def run_collaborative_server(host: str = "localhost", port: int = 8765,
                           enable_notifications: bool = True, enable_jupyter: bool = True):
    """Run the collaborative analysis server."""
    server = CollaborativeAnalysisServer(host, port, enable_notifications, enable_jupyter)

    try:
        asyncio.run(server.start_server())
    except KeyboardInterrupt:
        print("\n🛑 Server stopped")
    except Exception as e:
        print(f"❌ Server error: {e}")


def start_collaborative_session(project_name: str, config_file: str = "config/params.yaml") -> str:
    """Start a new collaborative session for a project."""
    server = CollaborativeAnalysisServer()
    session_id = server.create_session(project_name, config_file)

    print(f"🎯 Created collaborative session: {session_id}")
    print(f"🔗 Share this ID with collaborators: {session_id}")
    print(f"🛜 Connect to WebSocket: ws://{server.host}:{server.port}/{session_id}")

    # Create a sample notebook
    notebook_path = server.create_collaborative_notebook(session_id, "standard")
    if notebook_path:
        print(f"📓 Created collaborative notebook: {notebook_path}")

    return session_id


def join_collaborative_session(session_id: str, user_name: str, user_email: str) -> bool:
    """Join an existing collaborative session."""
    from scripts.collaboration_manager import User
    import uuid

    server = CollaborativeAnalysisServer()

    user = User(
        id=str(uuid.uuid4()),
        name=user_name,
        email=user_email,
        role="editor"
    )

    success = server.collaboration_manager.join_session(session_id, user)

    if success:
        print(f"✅ Joined session {session_id} as {user_name}")
        print(f"🛜 Connect to WebSocket: ws://{server.host}:{server.port}/{session_id}")

        # Get session info
        info = server.get_session_info(session_id)
        if info:
            print(f"📊 Session '{info['session']['name']}' - {len(info['session']['participants'])} participants")
    else:
        print(f"❌ Failed to join session {session_id}")

    return success


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNASEQ-MINI Collaborative Analysis Server")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Server command
    server_parser = subparsers.add_parser("server", help="Start collaborative server")
    server_parser.add_argument("--host", default="localhost", help="Server host")
    server_parser.add_argument("--port", type=int, default=8765, help="Server port")
    server_parser.add_argument("--no-notifications", action="store_true", help="Disable notifications")
    server_parser.add_argument("--no-jupyter", action="store_true", help="Disable Jupyter integration")

    # Create session command
    create_parser = subparsers.add_parser("create", help="Create a new collaborative session")
    create_parser.add_argument("project_name", help="Name of the analysis project")
    create_parser.add_argument("--config", default="config/params.yaml", help="Configuration file")

    # Join session command
    join_parser = subparsers.add_parser("join", help="Join an existing session")
    join_parser.add_argument("session_id", help="Session ID to join")
    join_parser.add_argument("user_name", help="Your name")
    join_parser.add_argument("user_email", help="Your email")

    # Setup notifications command
    setup_parser = subparsers.add_parser("setup-notifications", help="Setup notification channels")

    # List sessions command
    list_parser = subparsers.add_parser("list", help="List active sessions")

    # Session info command
    info_parser = subparsers.add_parser("info", help="Get session information")
    info_parser.add_argument("session_id", help="Session ID")

    args = parser.parse_args()

    if args.command == "server":
        run_collaborative_server(
            args.host, args.port,
            not args.no_notifications, not args.no_jupyter
        )

    elif args.command == "create":
        start_collaborative_session(args.project_name, args.config)

    elif args.command == "join":
        join_collaborative_session(args.session_id, args.user_name, args.user_email)

    elif args.command == "setup-notifications":
        server = CollaborativeAnalysisServer()
        server.setup_notifications()

    elif args.command == "list":
        server = CollaborativeAnalysisServer()
        sessions = server.list_active_sessions()
        if sessions:
            print("🗂️ Active Collaborative Sessions:")
            for session in sessions:
                print(f"  {session['id']}: {session['name']} ({session['participants']} participants)")
        else:
            print("📭 No active sessions")

    elif args.command == "info":
        server = CollaborativeAnalysisServer()
        info = server.get_session_info(args.session_id)
        if info:
            print(f"📊 Session: {info['session']['name']}")
            print(f"📍 Current step: {info['session']['current_step']}")
            print(f"👥 Participants: {len(info['session']['participants'])}")
            print(f"📅 Created: {info['session']['created_at']}")

            if info['recent_events']:
                print("
🕐 Recent Events:"                for event in info['recent_events'][:5]:
                    print(f"  {event['timestamp']}: {event['description']} ({event['user']})")
        else:
            print(f"❌ Session {args.session_id} not found")

    else:
        parser.print_help()










