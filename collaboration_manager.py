#!/usr/bin/env python3
"""
Real-Time Collaborative Analysis Environment for RNASEQ-MINI

This module enables multiple researchers to work together on RNA-seq analyses
with shared sessions, real-time updates, and collaborative workspaces.
"""

import json
import uuid
import time
import threading
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Any
import asyncio
import websockets
import sqlite3
from dataclasses import dataclass, asdict
import yaml
import os


@dataclass
class User:
    """Represents a user in the collaboration system."""
    id: str
    name: str
    email: str
    avatar_url: str = ""
    role: str = "viewer"  # viewer, editor, admin
    joined_at: datetime = None
    last_active: datetime = None

    def __post_init__(self):
        if self.joined_at is None:
            self.joined_at = datetime.now()
        if self.last_active is None:
            self.last_active = datetime.now()


@dataclass
class AnalysisSession:
    """Represents a shared analysis session."""
    id: str
    name: str
    description: str
    created_by: str
    created_at: datetime
    project_config: Dict[str, Any]
    status: str = "active"  # active, paused, completed, archived
    participants: Dict[str, User] = None
    current_step: str = "setup"
    shared_datasets: List[str] = None
    shared_results: List[str] = None

    def __post_init__(self):
        if self.participants is None:
            self.participants = {}
        if self.shared_datasets is None:
            self.shared_datasets = []
        if self.shared_results is None:
            self.shared_results = []


@dataclass
class AnalysisEvent:
    """Represents an event in the analysis workflow."""
    id: str
    session_id: str
    user_id: str
    event_type: str  # config_change, data_upload, analysis_start, result_view, etc.
    description: str
    data: Dict[str, Any]
    timestamp: datetime
    step: str = ""


class CollaborationDatabase:
    """SQLite database for storing collaboration data."""

    def __init__(self, db_path: str = ".collaboration.db"):
        self.db_path = Path(db_path)
        self.init_database()

    def init_database(self):
        """Initialize the collaboration database."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS sessions (
                    id TEXT PRIMARY KEY,
                    name TEXT NOT NULL,
                    description TEXT,
                    created_by TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    project_config TEXT,
                    status TEXT DEFAULT 'active',
                    current_step TEXT DEFAULT 'setup'
                )
            """)

            conn.execute("""
                CREATE TABLE IF NOT EXISTS participants (
                    session_id TEXT,
                    user_id TEXT,
                    name TEXT NOT NULL,
                    email TEXT,
                    avatar_url TEXT,
                    role TEXT DEFAULT 'viewer',
                    joined_at TEXT NOT NULL,
                    last_active TEXT NOT NULL,
                    PRIMARY KEY (session_id, user_id)
                )
            """)

            conn.execute("""
                CREATE TABLE IF NOT EXISTS events (
                    id TEXT PRIMARY KEY,
                    session_id TEXT NOT NULL,
                    user_id TEXT NOT NULL,
                    event_type TEXT NOT NULL,
                    description TEXT,
                    data TEXT,
                    timestamp TEXT NOT NULL,
                    step TEXT
                )
            """)

            conn.execute("""
                CREATE TABLE IF NOT EXISTS workspaces (
                    id TEXT PRIMARY KEY,
                    session_id TEXT NOT NULL,
                    name TEXT NOT NULL,
                    path TEXT NOT NULL,
                    type TEXT NOT NULL,  # dataset, config, results
                    created_by TEXT NOT NULL,
                    created_at TEXT NOT NULL
                )
            """)

    def create_session(self, session: AnalysisSession) -> str:
        """Create a new analysis session."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO sessions (id, name, description, created_by, created_at, project_config, status, current_step)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                session.id, session.name, session.description, session.created_by,
                session.created_at.isoformat(), json.dumps(session.project_config),
                session.status, session.current_step
            ))

            # Add creator as participant
            if session.participants:
                for user in session.participants.values():
                    conn.execute("""
                        INSERT OR REPLACE INTO participants
                        (session_id, user_id, name, email, avatar_url, role, joined_at, last_active)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        session.id, user.id, user.name, user.email, user.avatar_url,
                        user.role, user.joined_at.isoformat(), user.last_active.isoformat()
                    ))

        return session.id

    def get_session(self, session_id: str) -> Optional[AnalysisSession]:
        """Retrieve an analysis session."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT * FROM sessions WHERE id = ?", (session_id,))
            row = cursor.fetchone()

            if not row:
                return None

            # Get participants
            participants = {}
            cursor = conn.execute("SELECT * FROM participants WHERE session_id = ?", (session_id,))
            for p_row in cursor.fetchall():
                user = User(
                    id=p_row[1], name=p_row[2], email=p_row[3], avatar_url=p_row[4],
                    role=p_row[5], joined_at=datetime.fromisoformat(p_row[6]),
                    last_active=datetime.fromisoformat(p_row[7])
                )
                participants[user.id] = user

            return AnalysisSession(
                id=row[0], name=row[1], description=row[2], created_by=row[3],
                created_at=datetime.fromisoformat(row[4]), project_config=json.loads(row[5]),
                status=row[6], current_step=row[7], participants=participants
            )

    def add_participant(self, session_id: str, user: User):
        """Add a participant to a session."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT OR REPLACE INTO participants
                (session_id, user_id, name, email, avatar_url, role, joined_at, last_active)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                session_id, user.id, user.name, user.email, user.avatar_url,
                user.role, user.joined_at.isoformat(), user.last_active.isoformat()
            ))

    def log_event(self, event: AnalysisEvent):
        """Log an analysis event."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO events (id, session_id, user_id, event_type, description, data, timestamp, step)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                event.id, event.session_id, event.user_id, event.event_type,
                event.description, json.dumps(event.data), event.timestamp.isoformat(), event.step
            ))

    def get_session_events(self, session_id: str, limit: int = 100) -> List[AnalysisEvent]:
        """Get recent events for a session."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT * FROM events WHERE session_id = ? ORDER BY timestamp DESC LIMIT ?
            """, (session_id, limit))

            events = []
            for row in cursor.fetchall():
                event = AnalysisEvent(
                    id=row[0], session_id=row[1], user_id=row[2], event_type=row[3],
                    description=row[4], data=json.loads(row[5]), timestamp=datetime.fromisoformat(row[6]),
                    step=row[7]
                )
                events.append(event)

            return events


class CollaborationManager:
    """Main collaboration manager for real-time analysis."""

    def __init__(self, db_path: str = ".collaboration.db"):
        self.db = CollaborationDatabase(db_path)
        self.active_sessions: Dict[str, AnalysisSession] = {}
        self.websocket_connections: Dict[str, List[websockets.WebSocketServerProtocol]] = {}
        self.notification_queue: asyncio.Queue = asyncio.Queue()
        self.notification_manager = None  # Will be set by integration

        # Start background tasks
        self._start_background_tasks()

    def _start_background_tasks(self):
        """Start background maintenance tasks."""
        # Clean up inactive sessions periodically
        def cleanup_task():
            while True:
                time.sleep(300)  # Run every 5 minutes
                self._cleanup_inactive_sessions()

        cleanup_thread = threading.Thread(target=cleanup_task, daemon=True)
        cleanup_thread.start()

    def _cleanup_inactive_sessions(self):
        """Clean up sessions that have been inactive for too long."""
        cutoff = datetime.now() - timedelta(hours=24)

        with sqlite3.connect(self.db.db_path) as conn:
            cursor = conn.execute("""
                SELECT session_id FROM participants
                WHERE last_active < ?
                GROUP BY session_id
                HAVING COUNT(*) = (SELECT COUNT(*) FROM participants p2 WHERE p2.session_id = participants.session_id)
            """, (cutoff.isoformat(),))

            inactive_sessions = [row[0] for row in cursor.fetchall()]

            for session_id in inactive_sessions:
                conn.execute("UPDATE sessions SET status = 'inactive' WHERE id = ?", (session_id,))

    def create_session(self, name: str, description: str, created_by: str,
                      project_config: Dict[str, Any]) -> str:
        """Create a new collaborative analysis session."""
        session = AnalysisSession(
            id=str(uuid.uuid4()),
            name=name,
            description=description,
            created_by=created_by,
            created_at=datetime.now(),
            project_config=project_config
        )

        session_id = self.db.create_session(session)
        self.active_sessions[session_id] = session

        # Log creation event
        event = AnalysisEvent(
            id=str(uuid.uuid4()),
            session_id=session_id,
            user_id=created_by,
            event_type="session_created",
            description=f"Created analysis session: {name}",
            data={"name": name, "description": description},
            timestamp=datetime.now(),
            step="setup"
        )
        self.db.log_event(event)

        # Notify participants
        asyncio.create_task(self._broadcast_event(session_id, event))

        return session_id

    def join_session(self, session_id: str, user: User) -> bool:
        """Add a user to an existing session."""
        session = self.db.get_session(session_id)
        if not session or session.status not in ["active"]:
            return False

        session.participants[user.id] = user
        self.db.add_participant(session_id, user)
        self.active_sessions[session_id] = session

        # Log join event
        event = AnalysisEvent(
            id=str(uuid.uuid4()),
            session_id=session_id,
            user_id=user.id,
            event_type="user_joined",
            description=f"{user.name} joined the session",
            data={"user_name": user.name, "role": user.role},
            timestamp=datetime.now(),
            step=session.current_step
        )
        self.db.log_event(event)

        # Notify all participants
        asyncio.create_task(self._broadcast_event(session_id, event))

        return True

    def update_session_step(self, session_id: str, step: str, user_id: str):
        """Update the current analysis step."""
        session = self.active_sessions.get(session_id)
        if not session:
            return

        old_step = session.current_step
        session.current_step = step

        # Update in database
        with sqlite3.connect(self.db.db_path) as conn:
            conn.execute("UPDATE sessions SET current_step = ? WHERE id = ?", (step, session_id))

        # Log step change
        event = AnalysisEvent(
            id=str(uuid.uuid4()),
            session_id=session_id,
            user_id=user_id,
            event_type="step_changed",
            description=f"Analysis step changed: {old_step} → {step}",
            data={"old_step": old_step, "new_step": step},
            timestamp=datetime.now(),
            step=step
        )
        self.db.log_event(event)

        # Notify participants
        asyncio.create_task(self._broadcast_event(session_id, event))

    def log_analysis_event(self, session_id: str, user_id: str, event_type: str,
                          description: str, data: Dict[str, Any] = None, step: str = ""):
        """Log an analysis event."""
        if data is None:
            data = {}

        event = AnalysisEvent(
            id=str(uuid.uuid4()),
            session_id=session_id,
            user_id=user_id,
            event_type=event_type,
            description=description,
            data=data,
            timestamp=datetime.now(),
            step=step
        )

        self.db.log_event(event)

        # Notify participants via WebSocket
        asyncio.create_task(self._broadcast_event(session_id, event))

        # Send external notifications if manager is configured
        if self.notification_manager:
            asyncio.create_task(self._send_external_notification(event))

    def set_notification_manager(self, notification_manager):
        """Set the notification manager for external notifications."""
        self.notification_manager = notification_manager

    async def _send_external_notification(self, event: AnalysisEvent):
        """Send external notification for important events."""
        if not self.notification_manager:
            return

        # Format message based on event type
        message = self._format_event_message(event)

        # Send notification
        await self.notification_manager.send_notification(
            event.event_type,
            message,
            {
                "event_type": event.event_type,
                "session_id": event.session_id,
                "user_id": event.user_id,
                "step": event.step,
                "timestamp": event.timestamp.isoformat(),
                **event.data
            }
        )

    def _format_event_message(self, event: AnalysisEvent) -> str:
        """Format a user-friendly message for the event."""
        messages = {
            "session_created": f"🆕 New analysis session created: {event.data.get('name', 'Untitled')}",
            "user_joined": f"👋 {event.data.get('user_name', 'Someone')} joined the analysis session",
            "user_left": f"👋 {event.data.get('user_name', 'Someone')} left the analysis session",
            "step_changed": f"🔄 Analysis step changed: {event.data.get('old_step', 'unknown')} → {event.data.get('new_step', 'unknown')}",
            "analysis_started": "🚀 Analysis pipeline started",
            "analysis_completed": "✅ Analysis completed successfully",
            "analysis_failed": "❌ Analysis failed - check logs for details",
            "quality_gate_passed": "✅ Quality gate passed",
            "quality_gate_failed": "⚠️ Quality gate failed - review results",
            "data_uploaded": f"📁 New data uploaded: {event.data.get('file_name', 'unknown file')}",
            "config_changed": f"⚙️ Configuration updated: {event.description}",
            "custom_analysis": f"🔬 Custom analysis performed: {event.description}"
        }

        return messages.get(event.event_type, f"📢 {event.description}")

    async def _broadcast_event(self, session_id: str, event: AnalysisEvent):
        """Broadcast an event to all session participants."""
        if session_id not in self.websocket_connections:
            return

        event_data = {
            "type": "event",
            "event": asdict(event)
        }

        # Send to all connected websockets for this session
        disconnected = []
        for ws in self.websocket_connections[session_id]:
            try:
                await ws.send(json.dumps(event_data))
            except websockets.exceptions.ConnectionClosed:
                disconnected.append(ws)

        # Clean up disconnected websockets
        for ws in disconnected:
            self.websocket_connections[session_id].remove(ws)

    async def handle_websocket(self, websocket, path):
        """Handle WebSocket connections for real-time collaboration."""
        session_id = path.strip("/")

        if session_id not in self.websocket_connections:
            self.websocket_connections[session_id] = []

        self.websocket_connections[session_id].append(websocket)

        try:
            # Send current session state
            session = self.db.get_session(session_id)
            if session:
                await websocket.send(json.dumps({
                    "type": "session_state",
                    "session": asdict(session)
                }))

                # Send recent events
                events = self.db.get_session_events(session_id, 10)
                for event in events:
                    await websocket.send(json.dumps({
                        "type": "event",
                        "event": asdict(event)
                    }))

            # Keep connection alive and handle messages
            async for message in websocket:
                data = json.loads(message)

                if data.get("type") == "step_update":
                    self.update_session_step(
                        session_id,
                        data["step"],
                        data["user_id"]
                    )

                elif data.get("type") == "event_log":
                    self.log_analysis_event(
                        session_id,
                        data["user_id"],
                        data["event_type"],
                        data["description"],
                        data.get("data", {}),
                        data.get("step", "")
                    )

        except websockets.exceptions.ConnectionClosed:
            pass
        finally:
            if session_id in self.websocket_connections:
                if websocket in self.websocket_connections[session_id]:
                    self.websocket_connections[session_id].remove(websocket)

    def start_collaboration_server(self, host: str = "localhost", port: int = 8765):
        """Start the WebSocket server for real-time collaboration."""
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        start_server = websockets.serve(self.handle_websocket, host, port)

        print(f"🛜 Collaboration server starting on ws://{host}:{port}")
        print("Connect with: ws://localhost:8765/<session_id>")

        loop.run_until_complete(start_server)
        loop.run_forever()

    def export_session_data(self, session_id: str, output_dir: str = "exports"):
        """Export all session data for archival or sharing."""
        session = self.db.get_session(session_id)
        if not session:
            return False

        output_path = Path(output_dir) / f"session_{session_id}"
        output_path.mkdir(exist_ok=True)

        # Export session metadata
        session_data = {
            "session": asdict(session),
            "events": [asdict(event) for event in self.db.get_session_events(session_id)],
            "exported_at": datetime.now().isoformat()
        }

        with open(output_path / "session_data.json", 'w') as f:
            json.dump(session_data, f, indent=2, default=str)

        # Export configurations
        config_dir = output_path / "configs"
        config_dir.mkdir(exist_ok=True)

        # Copy shared configuration files
        for file_path in session.shared_datasets + session.shared_results:
            if Path(file_path).exists():
                shutil.copy2(file_path, config_dir)

        print(f"📤 Exported session data to {output_path}")
        return True


class WorkspaceManager:
    """Manages shared workspaces for collaborative analysis."""

    def __init__(self, base_dir: str = "shared_workspaces"):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(exist_ok=True)

    def create_workspace(self, session_id: str, name: str, created_by: str) -> str:
        """Create a shared workspace directory."""
        workspace_id = str(uuid.uuid4())
        workspace_path = self.base_dir / session_id / workspace_id
        workspace_path.mkdir(parents=True, exist_ok=True)

        # Create workspace metadata
        metadata = {
            "id": workspace_id,
            "session_id": session_id,
            "name": name,
            "created_by": created_by,
            "created_at": datetime.now().isoformat(),
            "path": str(workspace_path)
        }

        with open(workspace_path / ".workspace.json", 'w') as f:
            json.dump(metadata, f, indent=2)

        return workspace_id

    def share_file(self, workspace_id: str, file_path: str, file_type: str = "dataset"):
        """Share a file in a workspace."""
        workspace_path = self._get_workspace_path(workspace_id)
        if not workspace_path:
            return False

        file_name = Path(file_path).name
        shared_path = workspace_path / file_name

        # Copy file to workspace
        shutil.copy2(file_path, shared_path)

        # Update workspace metadata
        metadata = json.load(open(workspace_path / ".workspace.json"))
        if "shared_files" not in metadata:
            metadata["shared_files"] = []

        metadata["shared_files"].append({
            "name": file_name,
            "type": file_type,
            "shared_at": datetime.now().isoformat(),
            "original_path": file_path
        })

        with open(workspace_path / ".workspace.json", 'w') as f:
            json.dump(metadata, f, indent=2)

        return True

    def _get_workspace_path(self, workspace_id: str) -> Optional[Path]:
        """Get workspace path from ID."""
        for session_dir in self.base_dir.iterdir():
            if session_dir.is_dir():
                workspace_dir = session_dir / workspace_id
                if workspace_dir.exists():
                    return workspace_dir
        return None


def create_collaborative_session(project_name: str, config_file: str = "config/params.yaml") -> str:
    """Create a new collaborative analysis session."""
    # Load project configuration
    with open(config_file) as f:
        project_config = yaml.safe_load(f)

    manager = CollaborationManager()

    session_id = manager.create_session(
        name=f"Collaborative: {project_name}",
        description=f"Shared RNA-seq analysis for {project_name}",
        created_by="system",
        project_config=project_config
    )

    print(f"🎯 Created collaborative session: {session_id}")
    print(f"🔗 Share this ID with collaborators: {session_id}")
    print(f"🛜 Start server: python -m scripts.collaboration_manager server {session_id}")

    return session_id


def join_collaborative_session(session_id: str, user_name: str, user_email: str) -> bool:
    """Join an existing collaborative session."""
    manager = CollaborationManager()

    user = User(
        id=str(uuid.uuid4()),
        name=user_name,
        email=user_email,
        role="editor"
    )

    success = manager.join_session(session_id, user)

    if success:
        print(f"✅ Joined session {session_id} as {user_name}")
        print(f"🛜 Connect to WebSocket: ws://localhost:8765/{session_id}")
    else:
        print(f"❌ Failed to join session {session_id}")

    return success


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI Collaboration Manager")
    subparsers = parser.add_subparsers(dest="command")

    # Create session command
    create_parser = subparsers.add_parser("create", help="Create a new collaborative session")
    create_parser.add_argument("project_name", help="Name of the analysis project")
    create_parser.add_argument("--config", default="config/params.yaml", help="Configuration file")

    # Join session command
    join_parser = subparsers.add_parser("join", help="Join an existing session")
    join_parser.add_argument("session_id", help="Session ID to join")
    join_parser.add_argument("user_name", help="Your name")
    join_parser.add_argument("user_email", help="Your email")

    # Server command
    server_parser = subparsers.add_parser("server", help="Start collaboration server")
    server_parser.add_argument("--host", default="localhost", help="Server host")
    server_parser.add_argument("--port", type=int, default=8765, help="Server port")

    args = parser.parse_args()

    if args.command == "create":
        create_collaborative_session(args.project_name, args.config)
    elif args.command == "join":
        join_collaborative_session(args.session_id, args.user_name, args.user_email)
    elif args.command == "server":
        manager = CollaborationManager()
        manager.start_collaboration_server(args.host, args.port)
    else:
        parser.print_help()
