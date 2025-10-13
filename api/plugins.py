#!/usr/bin/env python3
"""
Plugin architecture for RNASEQ-MINI - enables extensibility and custom functionality.
Provides a framework for third-party modules and custom analysis tools.
"""

import json
import asyncio
import importlib
import importlib.util
import inspect
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Callable, Union
from datetime import datetime
import sys
import os

logger = logging.getLogger(__name__)


class Plugin:
    """Represents a plugin module."""

    def __init__(self, name: str, version: str, description: str,
                 author: str = "", license: str = "", dependencies: List[str] = None,
                 entry_points: Dict[str, str] = None):
        self.name = name
        self.version = version
        self.description = description
        self.author = author
        self.license = license
        self.dependencies = dependencies or []
        self.entry_points = entry_points or {}
        self.enabled = True
        self.loaded = False
        self.module = None
        self.created_at = datetime.now()

    def load(self) -> bool:
        """Load the plugin module."""
        try:
            # Import the plugin module
            if self.name in sys.modules:
                self.module = sys.modules[self.name]
            else:
                # Try to import from plugins directory
                plugin_path = Path("plugins") / self.name / "__init__.py"
                if plugin_path.exists():
                    spec = importlib.util.spec_from_file_location(self.name, plugin_path)
                    if spec and spec.loader:
                        self.module = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(self.module)
                        sys.modules[self.name] = self.module

            if self.module:
                # Check for required entry points
                required_methods = ['execute', 'validate']
                for method in required_methods:
                    if not hasattr(self.module, method):
                        logger.error(f"Plugin {self.name} missing required method: {method}")
                        return False

                self.loaded = True
                logger.info(f"Successfully loaded plugin: {self.name}")
                return True

        except Exception as e:
            logger.error(f"Error loading plugin {self.name}: {e}")
            return False

    def unload(self) -> bool:
        """Unload the plugin module."""
        try:
            if self.name in sys.modules:
                del sys.modules[self.name]

            self.module = None
            self.loaded = False
            logger.info(f"Unloaded plugin: {self.name}")
            return True

        except Exception as e:
            logger.error(f"Error unloading plugin {self.name}: {e}")
            return False

    async def execute(self, method_name: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a plugin method."""
        if not self.loaded or not self.module:
            raise RuntimeError(f"Plugin {self.name} is not loaded")

        if not hasattr(self.module, method_name):
            raise AttributeError(f"Plugin {self.name} does not have method: {method_name}")

        method = getattr(self.module, method_name)

        try:
            # Call the method (async or sync)
            if asyncio.iscoroutinefunction(method):
                result = await method(parameters)
            else:
                result = method(parameters)

            return {
                "success": True,
                "plugin": self.name,
                "method": method_name,
                "result": result,
                "timestamp": datetime.now().isoformat()
            }

        except Exception as e:
            logger.error(f"Error executing plugin method {self.name}.{method_name}: {e}")
            return {
                "success": False,
                "plugin": self.name,
                "method": method_name,
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            }


class PluginManager:
    """Manages plugin loading, execution, and lifecycle."""

    def __init__(self, plugins_dir: str = "plugins"):
        self.plugins_dir = Path(plugins_dir)
        self.plugins: Dict[str, Plugin] = {}
        self.event_handlers: Dict[str, List[Callable]] = {}

        # Create plugins directory if it doesn't exist
        self.plugins_dir.mkdir(exist_ok=True)

        # Load plugin registry
        self._load_plugin_registry()

    def _load_plugin_registry(self):
        """Load plugin registry from configuration."""
        registry_file = self.plugins_dir / "registry.json"

        if registry_file.exists():
            try:
                with open(registry_file, 'r') as f:
                    registry = json.load(f)

                for name, plugin_config in registry.items():
                    plugin = Plugin(**plugin_config)
                    self.plugins[name] = plugin

                logger.info(f"Loaded plugin registry with {len(self.plugins)} plugins")

            except Exception as e:
                logger.error(f"Error loading plugin registry: {e}")

    def _save_plugin_registry(self):
        """Save plugin registry to configuration."""
        try:
            registry = {}
            for name, plugin in self.plugins.items():
                registry[name] = {
                    "name": plugin.name,
                    "version": plugin.version,
                    "description": plugin.description,
                    "author": plugin.author,
                    "license": plugin.license,
                    "dependencies": plugin.dependencies,
                    "entry_points": plugin.entry_points,
                    "enabled": plugin.enabled,
                    "created_at": plugin.created_at.isoformat()
                }

            registry_file = self.plugins_dir / "registry.json"
            with open(registry_file, 'w') as f:
                json.dump(registry, f, indent=2)

        except Exception as e:
            logger.error(f"Error saving plugin registry: {e}")

    def register_plugin(self, plugin: Plugin) -> bool:
        """Register a new plugin."""
        try:
            self.plugins[plugin.name] = plugin
            self._save_plugin_registry()

            logger.info(f"Registered plugin: {plugin.name}")
            return True

        except Exception as e:
            logger.error(f"Error registering plugin {plugin.name}: {e}")
            return False

    def unregister_plugin(self, name: str) -> bool:
        """Unregister a plugin."""
        if name in self.plugins:
            plugin = self.plugins[name]

            # Unload if loaded
            if plugin.loaded:
                plugin.unload()

            # Remove from registry
            del self.plugins[name]
            self._save_plugin_registry()

            logger.info(f"Unregistered plugin: {name}")
            return True

        return False

    def list_plugins(self) -> Dict[str, Dict[str, Any]]:
        """List all registered plugins."""
        plugins_info = {}

        for name, plugin in self.plugins.items():
            plugins_info[name] = {
                "name": plugin.name,
                "version": plugin.version,
                "description": plugin.description,
                "author": plugin.author,
                "license": plugin.license,
                "enabled": plugin.enabled,
                "loaded": plugin.loaded,
                "dependencies": plugin.dependencies,
                "entry_points": list(plugin.entry_points.keys()),
                "created_at": plugin.created_at.isoformat()
            }

        return plugins_info

    def load_plugin(self, name: str) -> bool:
        """Load a specific plugin."""
        if name not in self.plugins:
            logger.error(f"Plugin {name} not found in registry")
            return False

        plugin = self.plugins[name]

        if plugin.loaded:
            logger.warning(f"Plugin {name} is already loaded")
            return True

        return plugin.load()

    def unload_plugin(self, name: str) -> bool:
        """Unload a specific plugin."""
        if name not in self.plugins:
            logger.error(f"Plugin {name} not found in registry")
            return False

        plugin = self.plugins[name]

        if not plugin.loaded:
            logger.warning(f"Plugin {name} is not loaded")
            return True

        return plugin.unload()

    def load_all_plugins(self) -> Dict[str, bool]:
        """Load all registered plugins."""
        results = {}

        for name, plugin in self.plugins.items():
            if plugin.enabled:
                results[name] = plugin.load()
            else:
                results[name] = True  # Already "loaded" (disabled)

        successful = sum(1 for result in results.values() if result)
        logger.info(f"Loaded {successful}/{len(results)} plugins")

        return results

    def unload_all_plugins(self) -> Dict[str, bool]:
        """Unload all loaded plugins."""
        results = {}

        for name, plugin in self.plugins.items():
            if plugin.loaded:
                results[name] = plugin.unload()
            else:
                results[name] = True

        return results

    async def execute_plugin(self, name: str, method: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a plugin method."""
        if name not in self.plugins:
            raise ValueError(f"Plugin {name} not found")

        plugin = self.plugins[name]

        if not plugin.enabled:
            raise ValueError(f"Plugin {name} is disabled")

        if not plugin.loaded:
            if not plugin.load():
                raise RuntimeError(f"Failed to load plugin {name}")

        return await plugin.execute(method, parameters)

    def enable_plugin(self, name: str) -> bool:
        """Enable a plugin."""
        if name in self.plugins:
            self.plugins[name].enabled = True
            self._save_plugin_registry()
            logger.info(f"Enabled plugin: {name}")
            return True
        return False

    def disable_plugin(self, name: str) -> bool:
        """Disable a plugin."""
        if name in self.plugins:
            plugin = self.plugins[name]

            # Unload if loaded
            if plugin.loaded:
                plugin.unload()

            plugin.enabled = False
            self._save_plugin_registry()
            logger.info(f"Disabled plugin: {name}")
            return True
        return False

    def get_plugin_info(self, name: str) -> Optional[Dict[str, Any]]:
        """Get detailed information about a plugin."""
        if name not in self.plugins:
            return None

        plugin = self.plugins[name]

        # Get available methods from loaded module
        methods = []
        if plugin.loaded and plugin.module:
            for attr_name in dir(plugin.module):
                attr = getattr(plugin.module, attr_name)
                if callable(attr) and not attr_name.startswith('_'):
                    methods.append({
                        "name": attr_name,
                        "signature": str(inspect.signature(attr)) if hasattr(attr, '__call__') else "Unknown"
                    })

        return {
            "name": plugin.name,
            "version": plugin.version,
            "description": plugin.description,
            "author": plugin.author,
            "license": plugin.license,
            "enabled": plugin.enabled,
            "loaded": plugin.loaded,
            "dependencies": plugin.dependencies,
            "entry_points": plugin.entry_points,
            "methods": methods,
            "created_at": plugin.created_at.isoformat()
        }

    def add_event_handler(self, event_type: str, handler: Callable):
        """Add a global event handler."""
        if event_type not in self.event_handlers:
            self.event_handlers[event_type] = []

        self.event_handlers[event_type].append(handler)

    async def emit_plugin_event(self, event_type: str, plugin_name: str, data: Dict[str, Any]):
        """Emit event to plugin event handlers."""
        if event_type in self.event_handlers:
            for handler in self.event_handlers[event_type]:
                try:
                    if asyncio.iscoroutinefunction(handler):
                        await handler(event_type, plugin_name, data)
                    else:
                        handler(event_type, plugin_name, data)
                except Exception as e:
                    logger.error(f"Error in plugin event handler for {event_type}: {e}")


# Plugin base class for developers
class RNASEQMiniPlugin:
    """Base class for RNASEQ-MINI plugins."""

    def __init__(self, name: str, version: str, description: str):
        self.name = name
        self.version = version
        self.description = description
        self.dependencies = []

    def validate(self, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate plugin parameters.

        Args:
            parameters: Plugin execution parameters

        Returns:
            Validation result with success status and any errors
        """
        # Default implementation - no validation
        return {"valid": True, "errors": []}

    def execute(self, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute the main plugin functionality.

        Args:
            parameters: Plugin execution parameters

        Returns:
            Execution results
        """
        raise NotImplementedError("Plugin must implement execute method")

    def cleanup(self) -> None:
        """Cleanup after plugin execution."""
        pass

    def get_metadata(self) -> Dict[str, Any]:
        """Get plugin metadata."""
        return {
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "dependencies": self.dependencies
        }


# Plugin registry for built-in plugins
class BuiltInPlugins:
    """Registry of built-in RNASEQ-MINI plugins."""

    @staticmethod
    def get_example_plugin() -> Plugin:
        """Get example plugin for demonstration."""
        return Plugin(
            name="example_plugin",
            version="1.0.0",
            description="Example plugin demonstrating plugin architecture",
            author="RNASEQ-MINI Team",
            license="MIT",
            dependencies=[],
            entry_points={
                "example_analysis": "example_plugin.ExamplePlugin.execute",
                "example_validation": "example_plugin.ExamplePlugin.validate"
            }
        )


# Plugin installation utilities
class PluginInstaller:
    """Helper for installing and managing plugins."""

    def __init__(self, plugins_dir: str = "plugins"):
        self.plugins_dir = Path(plugins_dir)

    def install_from_git(self, repo_url: str, branch: str = "main") -> bool:
        """Install plugin from Git repository."""
        try:
            import subprocess

            # Clone repository
            plugin_name = repo_url.split('/')[-1].replace('.git', '')
            plugin_path = self.plugins_dir / plugin_name

            if plugin_path.exists():
                logger.warning(f"Plugin directory {plugin_name} already exists")
                return False

            # Clone repository
            cmd = [
                "git", "clone",
                "--branch", branch,
                "--depth", "1",
                repo_url,
                str(plugin_path)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                logger.info(f"Successfully installed plugin from {repo_url}")
                return True
            else:
                logger.error(f"Failed to install plugin: {result.stderr}")
                return False

        except Exception as e:
            logger.error(f"Error installing plugin from git: {e}")
            return False

    def install_from_pip(self, package_name: str) -> bool:
        """Install plugin from PyPI."""
        try:
            import subprocess

            # Install package
            cmd = ["pip", "install", package_name]
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                logger.info(f"Successfully installed plugin package: {package_name}")
                return True
            else:
                logger.error(f"Failed to install plugin package: {result.stderr}")
                return False

        except Exception as e:
            logger.error(f"Error installing plugin from pip: {e}")
            return False

    def create_plugin_template(self, name: str, author: str = "", license: str = "MIT") -> bool:
        """Create a new plugin template."""
        try:
            plugin_dir = self.plugins_dir / name

            if plugin_dir.exists():
                logger.error(f"Plugin directory {name} already exists")
                return False

            # Create directory structure
            plugin_dir.mkdir(parents=True)

            # Create __init__.py
            init_content = f'''"""
{name} - Custom RNASEQ-MINI plugin.
"""

from .plugin import {name.title()}Plugin

__version__ = "0.1.0"
__author__ = "{author}"
__license__ = "{license}"
'''

            with open(plugin_dir / "__init__.py", 'w') as f:
                f.write(init_content)

            # Create plugin.py
            plugin_content = f'''"""
{name} plugin implementation.
"""

from api.plugins import RNASEQMiniPlugin


class {name.title()}Plugin(RNASEQMiniPlugin):
    """Custom plugin for {name}."""

    def __init__(self):
        super().__init__(
            name="{name}",
            version="0.1.0",
            description="Custom plugin description"
        )

    def validate(self, parameters):
        """Validate plugin parameters."""
        return {{"valid": True, "errors": []}}

    def execute(self, parameters):
        """Execute plugin functionality."""
        # Implement your plugin logic here
        return {{
            "status": "completed",
            "result": "Plugin executed successfully",
            "parameters": parameters
        }}

    def cleanup(self):
        """Cleanup after execution."""
        pass
'''

            with open(plugin_dir / "plugin.py", 'w') as f:
                f.write(plugin_content)

            # Create README.md
            readme_content = f'''# {name.title()} Plugin

Custom RNASEQ-MINI plugin for {name}.

## Installation

1. Place this directory in the `plugins/` folder of your RNASEQ-MINI installation
2. Register the plugin using the API or configuration

## Usage

```python
from api.plugins import PluginManager

manager = PluginManager()
manager.load_plugin("{name}")
result = await manager.execute_plugin("{name}", "execute", {{"param1": "value1"}})
```

## Configuration

Add to your `config/params.yaml`:

```yaml
plugins:
  {name}:
    enabled: true
    config:
      # Plugin-specific configuration
```
'''

            with open(plugin_dir / "README.md", 'w') as f:
                f.write(readme_content)

            logger.info(f"Created plugin template: {name}")
            return True

        except Exception as e:
            logger.error(f"Error creating plugin template: {e}")
            return False


# Plugin configuration management
class PluginConfig:
    """Plugin configuration manager."""

    def __init__(self, config_dir: str = "config"):
        self.config_dir = Path(config_dir)
        self.plugin_configs = {}

        # Load plugin configurations
        self._load_plugin_configs()

    def _load_plugin_configs(self):
        """Load plugin configurations."""
        config_file = self.config_dir / "plugins.yaml"

        if config_file.exists():
            try:
                import yaml
                with open(config_file, 'r') as f:
                    self.plugin_configs = yaml.safe_load(f) or {}

            except Exception as e:
                logger.error(f"Error loading plugin configurations: {e}")

    def get_plugin_config(self, plugin_name: str) -> Dict[str, Any]:
        """Get configuration for a specific plugin."""
        return self.plugin_configs.get(plugin_name, {})

    def set_plugin_config(self, plugin_name: str, config: Dict[str, Any]) -> bool:
        """Set configuration for a plugin."""
        try:
            self.plugin_configs[plugin_name] = config

            # Save to file
            config_file = self.config_dir / "plugins.yaml"
            import yaml

            with open(config_file, 'w') as f:
                yaml.dump(self.plugin_configs, f, default_flow_style=False)

            return True

        except Exception as e:
            logger.error(f"Error setting plugin config for {plugin_name}: {e}")
            return False

    def list_plugin_configs(self) -> Dict[str, Dict[str, Any]]:
        """List all plugin configurations."""
        return self.plugin_configs.copy()
