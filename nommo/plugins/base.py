"""
Base plugin interface for Nommo Engine.

Provides extension points for:
- Custom force fields
- Analysis algorithms
- Visualization backends
- Export formats
- Reaction rules
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import numpy as np

from nommo.utils.logging import get_logger

logger = get_logger("plugins")


class Plugin(ABC):
    """Base class for all plugins."""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Unique plugin name."""
        pass
    
    @property
    @abstractmethod
    def version(self) -> str:
        """Plugin version."""
        pass
    
    @property
    def description(self) -> str:
        """Plugin description."""
        return ""
    
    def initialize(self, config: Dict[str, Any]) -> None:
        """Initialize plugin with configuration."""
        pass
    
    def shutdown(self) -> None:
        """Clean up plugin resources."""
        pass


class ForceFieldPlugin(Plugin):
    """Base class for custom force field implementations."""
    
    @abstractmethod
    def calculate_forces(
        self,
        positions: np.ndarray,
        types: np.ndarray,
        box_size: np.ndarray,
        neighbor_list: Optional[List[Tuple[int, int]]] = None
    ) -> Tuple[np.ndarray, float]:
        """
        Calculate forces on all particles.
        
        Args:
            positions: Particle positions (N x 3)
            types: Particle type indices (N,)
            box_size: Simulation box dimensions
            neighbor_list: Optional pre-computed neighbor pairs
            
        Returns:
            forces: Force vectors (N x 3)
            potential_energy: Total potential energy
        """
        pass
    
    @abstractmethod
    def get_parameters(self) -> Dict[str, Any]:
        """Get force field parameters."""
        pass
    
    @abstractmethod
    def set_parameters(self, params: Dict[str, Any]) -> None:
        """Set force field parameters."""
        pass


class AnalysisPlugin(Plugin):
    """Base class for analysis algorithms."""
    
    @abstractmethod
    def analyze(
        self,
        universe: Any,
        history: Optional[List[Any]] = None
    ) -> Dict[str, Any]:
        """
        Perform analysis on universe state.
        
        Args:
            universe: Current universe state
            history: Optional historical states
            
        Returns:
            Analysis results dictionary
        """
        pass
    
    def supports_streaming(self) -> bool:
        """Whether this analysis supports streaming/online updates."""
        return False
    
    def update(self, universe: Any) -> Optional[Dict[str, Any]]:
        """Update analysis with new data (for streaming analysis)."""
        if not self.supports_streaming():
            raise NotImplementedError(f"{self.name} does not support streaming")
        return None


class VisualizationPlugin(Plugin):
    """Base class for visualization backends."""
    
    @abstractmethod
    def render(
        self,
        positions: np.ndarray,
        types: np.ndarray,
        bonds: Optional[List[Tuple[int, int]]] = None,
        colors: Optional[np.ndarray] = None,
        **kwargs
    ) -> Any:
        """
        Render particle system.
        
        Args:
            positions: Particle positions (N x 3)
            types: Particle types (N,)
            bonds: Optional bond pairs
            colors: Optional particle colors
            **kwargs: Additional rendering options
            
        Returns:
            Rendered output (format depends on backend)
        """
        pass
    
    def supports_animation(self) -> bool:
        """Whether this backend supports animation."""
        return False
    
    def create_animation(
        self,
        trajectory: List[np.ndarray],
        **kwargs
    ) -> Any:
        """Create animation from trajectory."""
        if not self.supports_animation():
            raise NotImplementedError(f"{self.name} does not support animation")
        return None


class ExportPlugin(Plugin):
    """Base class for data export formats."""
    
    @abstractmethod
    def export(
        self,
        data: Dict[str, Any],
        output_path: Path,
        **kwargs
    ) -> None:
        """
        Export data to file.
        
        Args:
            data: Data to export
            output_path: Output file path
            **kwargs: Format-specific options
        """
        pass
    
    @abstractmethod
    def import_data(
        self,
        input_path: Path,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Import data from file.
        
        Args:
            input_path: Input file path
            **kwargs: Format-specific options
            
        Returns:
            Imported data dictionary
        """
        pass
    
    @property
    @abstractmethod
    def file_extension(self) -> str:
        """File extension for this format."""
        pass


class ReactionPlugin(Plugin):
    """Base class for custom reaction rules."""
    
    @abstractmethod
    def can_react(
        self,
        particle1: Any,
        particle2: Any,
        distance: float,
        temperature: float
    ) -> bool:
        """
        Check if two particles can react.
        
        Args:
            particle1: First particle
            particle2: Second particle
            distance: Distance between particles
            temperature: System temperature
            
        Returns:
            Whether reaction can occur
        """
        pass
    
    @abstractmethod
    def apply_reaction(
        self,
        particle1: Any,
        particle2: Any,
        universe: Any
    ) -> List[Any]:
        """
        Apply reaction between particles.
        
        Args:
            particle1: First particle
            particle2: Second particle
            universe: Universe context
            
        Returns:
            List of product particles
        """
        pass
    
    @property
    @abstractmethod
    def reaction_type(self) -> str:
        """Type of reaction (e.g., 'synthesis', 'decomposition')."""
        pass


class PluginRegistry:
    """Registry for managing plugins."""
    
    def __init__(self):
        self._plugins: Dict[str, Plugin] = {}
        self._force_fields: Dict[str, ForceFieldPlugin] = {}
        self._analyzers: Dict[str, AnalysisPlugin] = {}
        self._visualizers: Dict[str, VisualizationPlugin] = {}
        self._exporters: Dict[str, ExportPlugin] = {}
        self._reactions: Dict[str, ReactionPlugin] = {}
        
    def register(self, plugin: Plugin) -> None:
        """Register a plugin."""
        if plugin.name in self._plugins:
            raise ValueError(f"Plugin {plugin.name} already registered")
            
        self._plugins[plugin.name] = plugin
        
        if isinstance(plugin, ForceFieldPlugin):
            self._force_fields[plugin.name] = plugin
        elif isinstance(plugin, AnalysisPlugin):
            self._analyzers[plugin.name] = plugin
        elif isinstance(plugin, VisualizationPlugin):
            self._visualizers[plugin.name] = plugin
        elif isinstance(plugin, ExportPlugin):
            self._exporters[plugin.name] = plugin
        elif isinstance(plugin, ReactionPlugin):
            self._reactions[plugin.name] = plugin
            
        logger.info(f"Registered plugin: {plugin.name} v{plugin.version}")
        
    def unregister(self, name: str) -> None:
        """Unregister a plugin."""
        if name not in self._plugins:
            raise ValueError(f"Plugin {name} not found")
            
        plugin = self._plugins.pop(name)
        
        self._force_fields.pop(name, None)
        self._analyzers.pop(name, None)
        self._visualizers.pop(name, None)
        self._exporters.pop(name, None)
        self._reactions.pop(name, None)
        
        plugin.shutdown()
        logger.info(f"Unregistered plugin: {name}")
        
    def get(self, name: str) -> Optional[Plugin]:
        """Get plugin by name."""
        return self._plugins.get(name)
    
    def get_force_field(self, name: str) -> Optional[ForceFieldPlugin]:
        """Get force field plugin."""
        return self._force_fields.get(name)
    
    def get_analyzer(self, name: str) -> Optional[AnalysisPlugin]:
        """Get analysis plugin."""
        return self._analyzers.get(name)
    
    def get_visualizer(self, name: str) -> Optional[VisualizationPlugin]:
        """Get visualization plugin."""
        return self._visualizers.get(name)
    
    def get_exporter(self, name: str) -> Optional[ExportPlugin]:
        """Get export plugin."""
        return self._exporters.get(name)
    
    def get_reaction(self, name: str) -> Optional[ReactionPlugin]:
        """Get reaction plugin."""
        return self._reactions.get(name)
    
    def list_plugins(self, plugin_type: Optional[str] = None) -> List[str]:
        """List registered plugins."""
        if plugin_type == "force_field":
            return list(self._force_fields.keys())
        elif plugin_type == "analysis":
            return list(self._analyzers.keys())
        elif plugin_type == "visualization":
            return list(self._visualizers.keys())
        elif plugin_type == "export":
            return list(self._exporters.keys())
        elif plugin_type == "reaction":
            return list(self._reactions.keys())
        else:
            return list(self._plugins.keys())
    
    def load_from_directory(self, directory: Path) -> None:
        """Load plugins from a directory."""
        import importlib.util
        import sys
        
        directory = Path(directory)
        if not directory.exists():
            logger.warning(f"Plugin directory not found: {directory}")
            return
            
        for path in directory.glob("*.py"):
            if path.stem.startswith("_"):
                continue
                
            spec = importlib.util.spec_from_file_location(path.stem, path)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                sys.modules[path.stem] = module
                spec.loader.exec_module(module)
                
                for attr_name in dir(module):
                    attr = getattr(module, attr_name)
                    if (
                        isinstance(attr, type) and
                        issubclass(attr, Plugin) and
                        attr != Plugin and
                        not attr.__name__.startswith("Base")
                    ):
                        try:
                            plugin = attr()
                            self.register(plugin)
                        except Exception as e:
                            logger.error(f"Failed to load plugin {attr.__name__}: {e}")


_registry = PluginRegistry()


def get_plugin_registry() -> PluginRegistry:
    """Get the global plugin registry."""
    return _registry