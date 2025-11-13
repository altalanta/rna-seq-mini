#!/usr/bin/env python3
"""
Real-Time Notification System for RNASEQ-MINI Collaborative Analysis

This module handles real-time notifications across multiple platforms
including Slack, Discord, email, and webhooks for collaborative sessions.
"""

import json
import asyncio
import aiohttp
import smtplib
from email.mime.text import MimeText
from email.mime.multipart import MimeMultipart
from typing import Dict, List, Optional, Any
import os
from datetime import datetime
from pathlib import Path
import yaml


class NotificationChannel:
    """Base class for notification channels."""

    def __init__(self, name: str, enabled: bool = True):
        self.name = name
        self.enabled = enabled

    async def send_notification(self, message: str, data: Dict[str, Any] = None):
        """Send a notification through this channel."""
        raise NotImplementedError


class SlackChannel(NotificationChannel):
    """Slack notification channel."""

    def __init__(self, webhook_url: str, channel: str = None, username: str = "RNASEQ-MINI"):
        super().__init__("slack")
        self.webhook_url = webhook_url
        self.channel = channel
        self.username = username

    async def send_notification(self, message: str, data: Dict[str, Any] = None):
        """Send notification to Slack."""
        if not self.enabled or not self.webhook_url:
            return

        payload = {
            "text": message,
            "username": self.username,
            "channel": self.channel
        }

        # Add attachments if data provided
        if data:
            attachment = {
                "color": self._get_color_for_event(data),
                "fields": []
            }

            for key, value in data.items():
                attachment["fields"].append({
                    "title": key.replace("_", " ").title(),
                    "value": str(value),
                    "short": True
                })

            payload["attachments"] = [attachment]

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(self.webhook_url, json=payload) as response:
                    if response.status != 200:
                        print(f"❌ Slack notification failed: {response.status}")
        except Exception as e:
            print(f"❌ Slack notification error: {e}")

    def _get_color_for_event(self, data: Dict[str, Any]) -> str:
        """Get color based on event type."""
        event_type = data.get("event_type", "")

        if event_type == "analysis_completed":
            return "good"
        elif event_type == "analysis_failed":
            return "danger"
        elif event_type == "quality_gate_failed":
            return "warning"
        elif event_type == "user_joined":
            return "good"
        else:
            return "#439FE0"


class DiscordChannel(NotificationChannel):
    """Discord notification channel."""

    def __init__(self, webhook_url: str, username: str = "RNASEQ-MINI"):
        super().__init__("discord")
        self.webhook_url = webhook_url
        self.username = username

    async def send_notification(self, message: str, data: Dict[str, Any] = None):
        """Send notification to Discord."""
        if not self.enabled or not self.webhook_url:
            return

        embed = {
            "title": "🔬 RNASEQ-MINI Update",
            "description": message,
            "timestamp": datetime.now().isoformat(),
            "color": self._get_color_for_event(data) if data else 0x439FE0
        }

        if data:
            embed["fields"] = []
            for key, value in data.items():
                embed["fields"].append({
                    "name": key.replace("_", " ").title(),
                    "value": str(value),
                    "inline": True
                })

        payload = {
            "username": self.username,
            "embeds": [embed]
        }

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(self.webhook_url, json=payload) as response:
                    if response.status != 204:
                        print(f"❌ Discord notification failed: {response.status}")
        except Exception as e:
            print(f"❌ Discord notification error: {e}")

    def _get_color_for_event(self, data: Dict[str, Any]) -> int:
        """Get color based on event type."""
        event_type = data.get("event_type", "")

        if event_type == "analysis_completed":
            return 0x00FF00  # Green
        elif event_type == "analysis_failed":
            return 0xFF0000  # Red
        elif event_type == "quality_gate_failed":
            return 0xFFFF00  # Yellow
        elif event_type == "user_joined":
            return 0x00FF00  # Green
        else:
            return 0x439FE0  # Blue


class EmailChannel(NotificationChannel):
    """Email notification channel."""

    def __init__(self, smtp_server: str, smtp_port: int, username: str, password: str,
                 from_email: str, to_emails: List[str]):
        super().__init__("email")
        self.smtp_server = smtp_server
        self.smtp_port = smtp_port
        self.username = username
        self.password = password
        self.from_email = from_email
        self.to_emails = to_emails

    async def send_notification(self, message: str, data: Dict[str, Any] = None):
        """Send notification via email."""
        if not self.enabled:
            return

        msg = MimeMultipart()
        msg['From'] = self.from_email
        msg['To'] = ", ".join(self.to_emails)
        msg['Subject'] = "🔬 RNASEQ-MINI Analysis Update"

        # Create HTML email
        html_content = f"""
        <html>
        <body>
            <h2>🔬 RNASEQ-MINI Analysis Update</h2>
            <p>{message}</p>
        """

        if data:
            html_content += """
            <h3>Details:</h3>
            <ul>
            """
            for key, value in data.items():
                html_content += f"<li><strong>{key.replace('_', ' ').title()}:</strong> {value}</li>\n"
            html_content += "</ul>"

        html_content += """
            <p>
                <small>
                    This notification was sent by RNASEQ-MINI collaboration system.<br>
                    To manage notifications, update your collaboration settings.
                </small>
            </p>
        </body>
        </html>
        """

        msg.attach(MimeText(html_content, 'html'))

        try:
            with smtplib.SMTP(self.smtp_server, self.smtp_port) as server:
                server.starttls()
                server.login(self.username, self.password)
                server.send_message(msg)
        except Exception as e:
            print(f"❌ Email notification error: {e}")


class WebhookChannel(NotificationChannel):
    """Generic webhook notification channel."""

    def __init__(self, webhook_url: str, headers: Dict[str, str] = None):
        super().__init__("webhook")
        self.webhook_url = webhook_url
        self.headers = headers or {"Content-Type": "application/json"}

    async def send_notification(self, message: str, data: Dict[str, Any] = None):
        """Send notification via webhook."""
        if not self.enabled or not self.webhook_url:
            return

        payload = {
            "message": message,
            "timestamp": datetime.now().isoformat(),
            "source": "rnaseq-mini"
        }

        if data:
            payload["data"] = data

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    self.webhook_url,
                    json=payload,
                    headers=self.headers
                ) as response:
                    if response.status not in [200, 201, 204]:
                        print(f"❌ Webhook notification failed: {response.status}")
        except Exception as e:
            print(f"❌ Webhook notification error: {e}")


class NotificationManager:
    """Main notification manager for the collaboration system."""

    def __init__(self, config_file: str = "config/notifications.yaml"):
        self.config_file = Path(config_file)
        self.channels: Dict[str, NotificationChannel] = {}
        self.event_filters: Dict[str, List[str]] = {}

        self.load_config()

    def load_config(self):
        """Load notification configuration."""
        if not self.config_file.exists():
            self._create_default_config()
            return

        with open(self.config_file) as f:
            config = yaml.safe_load(f)

        # Load notification channels
        for channel_name, channel_config in config.get("channels", {}).items():
            if not channel_config.get("enabled", True):
                continue

            channel_type = channel_config.get("type", "webhook")

            if channel_type == "slack":
                self.channels[channel_name] = SlackChannel(
                    webhook_url=channel_config["webhook_url"],
                    channel=channel_config.get("channel"),
                    username=channel_config.get("username", "RNASEQ-MINI")
                )
            elif channel_type == "discord":
                self.channels[channel_name] = DiscordChannel(
                    webhook_url=channel_config["webhook_url"],
                    username=channel_config.get("username", "RNASEQ-MINI")
                )
            elif channel_type == "email":
                self.channels[channel_name] = EmailChannel(
                    smtp_server=channel_config["smtp_server"],
                    smtp_port=channel_config["smtp_port"],
                    username=channel_config["username"],
                    password=channel_config["password"],
                    from_email=channel_config["from_email"],
                    to_emails=channel_config["to_emails"]
                )
            elif channel_type == "webhook":
                self.channels[channel_name] = WebhookChannel(
                    webhook_url=channel_config["webhook_url"],
                    headers=channel_config.get("headers", {})
                )

        # Load event filters
        self.event_filters = config.get("event_filters", {})

    def _create_default_config(self):
        """Create default notification configuration."""
        default_config = {
            "channels": {
                "default_webhook": {
                    "type": "webhook",
                    "enabled": False,
                    "webhook_url": "https://your-webhook-url.com/endpoint"
                }
            },
            "event_filters": {
                "*": ["analysis_completed", "analysis_failed", "quality_gate_failed", "user_joined", "user_left"]
            }
        }

        self.config_file.parent.mkdir(exist_ok=True)
        with open(self.config_file, 'w') as f:
            yaml.dump(default_config, f, default_flow_style=False, indent=2)

        print(f"📋 Created default notification config: {self.config_file}")
        print("💡 Edit this file to configure your notification channels")

    def add_channel(self, name: str, channel: NotificationChannel):
        """Add a notification channel."""
        self.channels[name] = channel
        self.save_config()

    def remove_channel(self, name: str):
        """Remove a notification channel."""
        if name in self.channels:
            del self.channels[name]
            self.save_config()

    def save_config(self):
        """Save current configuration to file."""
        config = {
            "channels": {},
            "event_filters": self.event_filters
        }

        for name, channel in self.channels.items():
            channel_config = {
                "enabled": channel.enabled,
                "type": channel.name
            }

            if isinstance(channel, SlackChannel):
                channel_config.update({
                    "webhook_url": channel.webhook_url,
                    "channel": channel.channel,
                    "username": channel.username
                })
            elif isinstance(channel, DiscordChannel):
                channel_config.update({
                    "webhook_url": channel.webhook_url,
                    "username": channel.username
                })
            elif isinstance(channel, EmailChannel):
                channel_config.update({
                    "smtp_server": channel.smtp_server,
                    "smtp_port": channel.smtp_port,
                    "username": channel.username,
                    "password": "***",  # Don't save actual password
                    "from_email": channel.from_email,
                    "to_emails": channel.to_emails
                })
            elif isinstance(channel, WebhookChannel):
                channel_config.update({
                    "webhook_url": channel.webhook_url,
                    "headers": channel.headers
                })

            config["channels"][name] = channel_config

        with open(self.config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=2)

    async def send_notification(self, event_type: str, message: str,
                               data: Dict[str, Any] = None, session_id: str = None):
        """Send notification to all enabled channels."""
        # Check if this event type should be filtered
        enabled_channels = []

        for channel_name, channel in self.channels.items():
            if not channel.enabled:
                continue

            # Check if channel should receive this event type
            channel_filters = self.event_filters.get(channel_name, self.event_filters.get("*", []))
            if event_type not in channel_filters and "*" not in channel_filters:
                continue

            enabled_channels.append(channel)

        if not enabled_channels:
            return

        # Send to all enabled channels
        tasks = []
        for channel in enabled_channels:
            task = channel.send_notification(message, data)
            tasks.append(task)

        if tasks:
            await asyncio.gather(*tasks, return_exceptions=True)

    def setup_slack_integration(self, webhook_url: str, channel: str = None):
        """Set up Slack integration."""
        slack = SlackChannel(webhook_url, channel)
        self.add_channel("slack", slack)
        print("✅ Slack integration configured")

    def setup_discord_integration(self, webhook_url: str):
        """Set up Discord integration."""
        discord = DiscordChannel(webhook_url)
        self.add_channel("discord", discord)
        print("✅ Discord integration configured")

    def setup_email_integration(self, smtp_server: str, smtp_port: int,
                               username: str, password: str, from_email: str, to_emails: List[str]):
        """Set up email integration."""
        email = EmailChannel(smtp_server, smtp_port, username, password, from_email, to_emails)
        self.add_channel("email", email)
        print("✅ Email integration configured")

    def test_notifications(self):
        """Test all notification channels."""
        async def test():
            test_data = {
                "event_type": "test_notification",
                "session_id": "test-session",
                "timestamp": datetime.now().isoformat()
            }

            await self.send_notification(
                "test_notification",
                "🧪 This is a test notification from RNASEQ-MINI",
                test_data
            )

        asyncio.run(test())
        print("✅ Test notifications sent")


# Global notification manager instance
_notification_manager: Optional[NotificationManager] = None


def get_notification_manager() -> NotificationManager:
    """Get the global notification manager instance."""
    global _notification_manager
    if _notification_manager is None:
        _notification_manager = NotificationManager()
    return _notification_manager


def setup_notifications():
    """Interactive setup for notification channels."""
    manager = get_notification_manager()

    print("🔔 RNASEQ-MINI Notification Setup")
    print("=" * 50)

    while True:
        print("\nAvailable notification channels:")
        print("1. Slack")
        print("2. Discord")
        print("3. Email")
        print("4. Webhook")
        print("5. Test notifications")
        print("6. Done")

        choice = input("\nSelect option (1-6): ").strip()

        if choice == "1":
            webhook_url = input("Slack webhook URL: ").strip()
            channel = input("Slack channel (optional): ").strip() or None
            if webhook_url:
                manager.setup_slack_integration(webhook_url, channel)

        elif choice == "2":
            webhook_url = input("Discord webhook URL: ").strip()
            if webhook_url:
                manager.setup_discord_integration(webhook_url)

        elif choice == "3":
            print("Email configuration:")
            smtp_server = input("SMTP server: ").strip()
            smtp_port = int(input("SMTP port (587): ").strip() or "587")
            username = input("Username: ").strip()
            password = input("Password: ").strip()
            from_email = input("From email: ").strip()
            to_emails_input = input("To emails (comma-separated): ").strip()
            to_emails = [email.strip() for email in to_emails_input.split(",") if email.strip()]

            if all([smtp_server, username, password, from_email, to_emails]):
                manager.setup_email_integration(smtp_server, smtp_port, username, password, from_email, to_emails)

        elif choice == "4":
            webhook_url = input("Webhook URL: ").strip()
            if webhook_url:
                webhook = WebhookChannel(webhook_url)
                name = input("Channel name: ").strip() or "webhook"
                manager.add_channel(name, webhook)
                print(f"✅ Webhook '{name}' configured")

        elif choice == "5":
            print("🧪 Testing notifications...")
            manager.test_notifications()

        elif choice == "6":
            break

        else:
            print("❌ Invalid option")

    print("✅ Notification setup complete!")
    print(f"📋 Configuration saved to: {manager.config_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI Notification Manager")
    subparsers = parser.add_subparsers(dest="command")

    # Setup command
    setup_parser = subparsers.add_parser("setup", help="Interactive notification setup")

    # Test command
    test_parser = subparsers.add_parser("test", help="Test notification channels")

    # Send command
    send_parser = subparsers.add_parser("send", help="Send a custom notification")
    send_parser.add_argument("event_type", help="Event type")
    send_parser.add_argument("message", help="Notification message")
    send_parser.add_argument("--data", help="JSON data for notification")

    args = parser.parse_args()

    if args.command == "setup":
        setup_notifications()
    elif args.command == "test":
        manager = get_notification_manager()
        manager.test_notifications()
    elif args.command == "send":
        manager = get_notification_manager()
        data = json.loads(args.data) if args.data else {}

        async def send():
            await manager.send_notification(args.event_type, args.message, data)

        asyncio.run(send())
    else:
        parser.print_help()







