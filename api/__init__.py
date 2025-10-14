"""
RNASEQ-MINI API - Programmatic interface for all pipeline operations.
Provides REST API, SDK, webhooks, and plugin architecture.
"""

from .server import RNASEQMiniAPI
from .client import RNASEQMiniClient
from .webhooks import WebhookManager
from .plugins import PluginManager
from .auth import APIAuthenticator

__all__ = [
    'RNASEQMiniAPI',
    'RNASEQMiniClient',
    'WebhookManager',
    'PluginManager',
    'APIAuthenticator'
]

__version__ = "1.0.0"


