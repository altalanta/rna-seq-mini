#!/usr/bin/env python3
"""
Webhook system for RNASEQ-MINI - enables external integrations and notifications.
Provides event-driven notifications for pipeline events, errors, and completions.
"""

import json
import asyncio
import aiohttp
import logging
from typing import Dict, List, Optional, Any, Callable
from pathlib import Path
from datetime import datetime
import hmac
import hashlib
import base64

logger = logging.getLogger(__name__)


class Webhook:
    """Represents a webhook endpoint."""

    def __init__(self, url: str, secret: str = None, events: List[str] = None,
                 headers: Dict[str, str] = None, enabled: bool = True):
        self.url = url
        self.secret = secret  # For HMAC signature verification
        self.events = events or ["*"]  # Events to subscribe to
        self.headers = headers or {}
        self.enabled = enabled
        self.created_at = datetime.now()
        self.last_triggered = None
        self.failure_count = 0

    def should_trigger(self, event_type: str) -> bool:
        """Check if webhook should trigger for given event."""
        return self.enabled and (event_type in self.events or "*" in self.events)

    def generate_signature(self, payload: str) -> str:
        """Generate HMAC signature for payload."""
        if not self.secret:
            return ""

        message = payload.encode('utf-8')
        secret = self.secret.encode('utf-8')
        signature = hmac.new(secret, message, hashlib.sha256).hexdigest()
        return f"sha256={signature}"

    async def send(self, event_type: str, payload: Dict[str, Any]) -> bool:
        """Send webhook payload to endpoint."""
        if not self.should_trigger(event_type):
            return True  # Not subscribed to this event

        try:
            # Prepare payload
            webhook_payload = {
                "event": event_type,
                "timestamp": datetime.now().isoformat(),
                "data": payload,
                "source": "rnaseq-mini"
            }

            # Convert to JSON string for signature
            payload_str = json.dumps(webhook_payload, sort_keys=True)

            # Prepare headers
            headers = {
                "Content-Type": "application/json",
                "User-Agent": "RNASEQ-MINI-Webhook/1.0",
                **self.headers
            }

            # Add signature if secret is configured
            if self.secret:
                signature = self.generate_signature(payload_str)
                headers["X-Hub-Signature-256"] = signature

            # Send request
            timeout = aiohttp.ClientTimeout(total=30)
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.post(
                    self.url,
                    json=webhook_payload,
                    headers=headers
                ) as response:
                    success = response.status < 400

                    if success:
                        self.last_triggered = datetime.now()
                        self.failure_count = 0
                        logger.info(f"Webhook sent successfully to {self.url} for event {event_type}")
                    else:
                        self.failure_count += 1
                        logger.error(f"Webhook failed for {self.url} (status {response.status}): {await response.text()}")

                    return success

        except Exception as e:
            self.failure_count += 1
            logger.error(f"Error sending webhook to {self.url}: {e}")
            return False


class WebhookManager:
    """Manages webhook registrations and triggers."""

    def __init__(self, config_file: str = "config/webhooks.json"):
        self.config_file = Path(config_file)
        self.webhooks: Dict[str, Webhook] = {}
        self.event_handlers: Dict[str, List[Callable]] = {}
        self._load_config()

    def _load_config(self):
        """Load webhook configuration from file."""
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    config = json.load(f)

                for name, webhook_config in config.items():
                    webhook = Webhook(**webhook_config)
                    self.webhooks[name] = webhook

                logger.info(f"Loaded {len(self.webhooks)} webhooks from {self.config_file}")

            except Exception as e:
                logger.error(f"Error loading webhook config: {e}")

    def _save_config(self):
        """Save webhook configuration to file."""
        try:
            config = {}
            for name, webhook in self.webhooks.items():
                config[name] = {
                    "url": webhook.url,
                    "secret": webhook.secret,
                    "events": webhook.events,
                    "headers": webhook.headers,
                    "enabled": webhook.enabled,
                    "created_at": webhook.created_at.isoformat(),
                    "last_triggered": webhook.last_triggered.isoformat() if webhook.last_triggered else None,
                    "failure_count": webhook.failure_count
                }

            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=2)

        except Exception as e:
            logger.error(f"Error saving webhook config: {e}")

    def register_webhook(self, name: str, url: str, secret: str = None,
                        events: List[str] = None, headers: Dict[str, str] = None) -> bool:
        """
        Register a new webhook.

        Args:
            name: Unique webhook identifier
            url: Webhook endpoint URL
            secret: Optional secret for signature verification
            events: Events to subscribe to
            headers: Additional headers to send

        Returns:
            Success status
        """
        try:
            webhook = Webhook(url, secret, events, headers)
            self.webhooks[name] = webhook
            self._save_config()

            logger.info(f"Registered webhook '{name}' for events: {webhook.events}")
            return True

        except Exception as e:
            logger.error(f"Error registering webhook {name}: {e}")
            return False

    def unregister_webhook(self, name: str) -> bool:
        """Unregister a webhook."""
        if name in self.webhooks:
            del self.webhooks[name]
            self._save_config()
            logger.info(f"Unregistered webhook '{name}'")
            return True

        return False

    def list_webhooks(self) -> Dict[str, Dict[str, Any]]:
        """List all registered webhooks."""
        webhooks_info = {}
        for name, webhook in self.webhooks.items():
            webhooks_info[name] = {
                "url": webhook.url,
                "events": webhook.events,
                "enabled": webhook.enabled,
                "created_at": webhook.created_at.isoformat(),
                "last_triggered": webhook.last_triggered.isoformat() if webhook.last_triggered else None,
                "failure_count": webhook.failure_count
            }

        return webhooks_info

    def enable_webhook(self, name: str) -> bool:
        """Enable a webhook."""
        if name in self.webhooks:
            self.webhooks[name].enabled = True
            self._save_config()
            logger.info(f"Enabled webhook '{name}'")
            return True
        return False

    def disable_webhook(self, name: str) -> bool:
        """Disable a webhook."""
        if name in self.webhooks:
            self.webhooks[name].enabled = False
            self._save_config()
            logger.info(f"Disabled webhook '{name}'")
            return True
        return False

    async def trigger_webhook(self, event_type: str, payload: Dict[str, Any]) -> Dict[str, Any]:
        """
        Trigger webhooks for a specific event.

        Args:
            event_type: Type of event (job_completed, job_failed, etc.)
            payload: Event payload data

        Returns:
            Results of webhook triggers
        """
        results = {
            "event_type": event_type,
            "timestamp": datetime.now().isoformat(),
            "webhooks_triggered": 0,
            "webhooks_failed": 0,
            "webhooks_skipped": 0,
            "details": []
        }

        # Trigger all relevant webhooks concurrently
        tasks = []
        for name, webhook in self.webhooks.items():
            if webhook.should_trigger(event_type):
                task = webhook.send(event_type, payload)
                tasks.append((name, task))

        # Execute all webhooks
        if tasks:
            webhook_results = await asyncio.gather(*[task for _, task in tasks], return_exceptions=True)

            for i, (name, result) in enumerate(zip([name for name, _ in tasks], webhook_results)):
                webhook = self.webhooks[[name for name, _ in tasks][i]]

                if isinstance(result, Exception):
                    results["webhooks_failed"] += 1
                    results["details"].append({
                        "webhook": name,
                        "status": "error",
                        "error": str(result)
                    })
                elif result:
                    results["webhooks_triggered"] += 1
                    results["details"].append({
                        "webhook": name,
                        "status": "success"
                    })
                else:
                    results["webhooks_failed"] += 1
                    results["details"].append({
                        "webhook": name,
                        "status": "failed"
                    })
        else:
            results["webhooks_skipped"] = len(self.webhooks)

        logger.info(f"Triggered {results['webhooks_triggered']} webhooks for event {event_type}")
        return results

    def add_event_handler(self, event_type: str, handler: Callable):
        """Add a local event handler."""
        if event_type not in self.event_handlers:
            self.event_handlers[event_type] = []

        self.event_handlers[event_type].append(handler)
        logger.info(f"Added event handler for {event_type}")

    async def emit_event(self, event_type: str, payload: Dict[str, Any]):
        """Emit an event to both webhooks and local handlers."""
        # Trigger webhooks
        webhook_results = await self.trigger_webhook(event_type, payload)

        # Trigger local handlers
        if event_type in self.event_handlers:
            for handler in self.event_handlers[event_type]:
                try:
                    if asyncio.iscoroutinefunction(handler):
                        await handler(event_type, payload)
                    else:
                        handler(event_type, payload)
                except Exception as e:
                    logger.error(f"Error in event handler for {event_type}: {e}")

        return webhook_results

    def get_failed_webhooks(self) -> List[str]:
        """Get list of webhooks with recent failures."""
        failed = []
        for name, webhook in self.webhooks.items():
            if webhook.failure_count > 0:
                failed.append(name)

        return failed

    def reset_webhook_failures(self, name: str) -> bool:
        """Reset failure count for a webhook."""
        if name in self.webhooks:
            self.webhooks[name].failure_count = 0
            self._save_config()
            logger.info(f"Reset failure count for webhook '{name}'")
            return True
        return False


class WebhookEventTypes:
    """Standard webhook event types."""

    # Pipeline events
    PIPELINE_STARTED = "pipeline_started"
    PIPELINE_COMPLETED = "pipeline_completed"
    PIPELINE_FAILED = "pipeline_failed"

    # Job events
    JOB_QUEUED = "job_queued"
    JOB_STARTED = "job_started"
    JOB_COMPLETED = "job_completed"
    JOB_FAILED = "job_failed"

    # Quality events
    QUALITY_ASSESSMENT_COMPLETED = "quality_assessment_completed"
    QUALITY_GATE_PASSED = "quality_gate_passed"
    QUALITY_GATE_FAILED = "quality_gate_failed"

    # Configuration events
    CONFIG_OPTIMIZED = "config_optimized"
    CONFIG_VALIDATED = "config_validated"

    # Multi-omics events
    MULTIOMICS_INTEGRATION_COMPLETED = "multiomics_integration_completed"
    MULTIOMICS_NORMALIZATION_COMPLETED = "multiomics_normalization_completed"

    # Cache events
    CACHE_CLEANUP_COMPLETED = "cache_cleanup_completed"

    # Plugin events
    PLUGIN_EXECUTED = "plugin_executed"
    PLUGIN_FAILED = "plugin_failed"

    # System events
    SYSTEM_STARTUP = "system_startup"
    SYSTEM_SHUTDOWN = "system_shutdown"
    BACKUP_COMPLETED = "backup_completed"


# Utility functions for webhook management
def create_slack_webhook(url: str, channel: str = None) -> Dict[str, str]:
    """Create headers for Slack webhook."""
    headers = {"Content-Type": "application/json"}

    if channel:
        # Custom channel routing
        headers["X-Slack-Channel"] = channel

    return headers


def create_discord_webhook(url: str, username: str = "RNASEQ-MINI") -> Dict[str, str]:
    """Create headers for Discord webhook."""
    headers = {"Content-Type": "application/json"}

    # Discord expects specific payload format
    return headers


def create_teams_webhook(url: str, title: str = "RNASEQ-MINI Notification") -> Dict[str, str]:
    """Create headers for Microsoft Teams webhook."""
    headers = {"Content-Type": "application/json"}

    # Teams expects specific payload format
    return headers


def validate_webhook_signature(payload: str, signature: str, secret: str) -> bool:
    """Validate webhook signature."""
    if not signature or not secret:
        return False

    try:
        # Extract signature value
        if signature.startswith("sha256="):
            expected_signature = signature[7:]  # Remove "sha256=" prefix
        else:
            expected_signature = signature

        # Generate signature for comparison
        message = payload.encode('utf-8')
        secret_bytes = secret.encode('utf-8')
        actual_signature = hmac.new(secret_bytes, message, hashlib.sha256).hexdigest()

        return hmac.compare_digest(expected_signature, actual_signature)

    except Exception:
        return False


# Example webhook payload formatters
def format_slack_message(event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
    """Format message for Slack webhook."""
    color_map = {
        "job_completed": "good",
        "job_failed": "danger",
        "pipeline_completed": "good",
        "pipeline_failed": "danger",
        "quality_gate_passed": "good",
        "quality_gate_failed": "warning"
    }

    color = color_map.get(event_type, "good")

    return {
        "attachments": [{
            "color": color,
            "title": f"RNASEQ-MINI: {event_type.replace('_', ' ').title()}",
            "fields": [
                {"title": "Event", "value": event_type, "short": True},
                {"title": "Timestamp", "value": data.get("timestamp", "Unknown"), "short": True}
            ],
            "footer": "RNASEQ-MINI",
            "ts": int(datetime.now().timestamp())
        }]
    }


def format_discord_message(event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
    """Format message for Discord webhook."""
    embed_color = {
        "job_completed": 0x00ff00,
        "job_failed": 0xff0000,
        "pipeline_completed": 0x00ff00,
        "pipeline_failed": 0xff0000,
        "quality_gate_passed": 0x00ff00,
        "quality_gate_failed": 0xffff00
    }.get(event_type, 0x00ff00)

    return {
        "embeds": [{
            "title": f"RNASEQ-MINI: {event_type.replace('_', ' ').title()}",
            "color": embed_color,
            "fields": [
                {"name": "Event", "value": event_type, "inline": True},
                {"name": "Timestamp", "value": data.get("timestamp", "Unknown"), "inline": True}
            ],
            "footer": {"text": "RNASEQ-MINI"},
            "timestamp": datetime.now().isoformat()
        }]
    }


def format_teams_message(event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
    """Format message for Microsoft Teams webhook."""
    return {
        "@type": "MessageCard",
        "@context": "http://schema.org/extensions",
        "themeColor": "0076D7",
        "summary": f"RNASEQ-MINI: {event_type.replace('_', ' ').title()}",
        "sections": [{
            "activityTitle": f"RNASEQ-MINI: {event_type.replace('_', ' ').title()}",
            "activitySubtitle": data.get("timestamp", "Unknown"),
            "facts": [
                {"name": "Event", "value": event_type},
                {"name": "Status", "value": data.get("status", "Unknown")}
            ]
        }]
    }


# Webhook payload formatters registry
WEBHOOK_FORMATTERS = {
    "slack": format_slack_message,
    "discord": format_discord_message,
    "teams": format_teams_message
}


class WebhookFormatter:
    """Helper class for formatting webhook payloads."""

    @staticmethod
    def format_for_platform(event_type: str, data: Dict[str, Any], platform: str) -> Dict[str, Any]:
        """Format payload for specific platform."""
        formatter = WEBHOOK_FORMATTERS.get(platform.lower())
        if formatter:
            return formatter(event_type, data)

        # Default JSON format
        return {
            "event": event_type,
            "timestamp": datetime.now().isoformat(),
            "data": data,
            "source": "rnaseq-mini"
        }

    @staticmethod
    def detect_platform_from_url(url: str) -> str:
        """Detect platform from webhook URL."""
        url_lower = url.lower()

        if "slack" in url_lower or "hooks.slack" in url_lower:
            return "slack"
        elif "discord" in url_lower or "discordapp" in url_lower:
            return "discord"
        elif "teams" in url_lower or "office" in url_lower:
            return "teams"

        return "generic"



