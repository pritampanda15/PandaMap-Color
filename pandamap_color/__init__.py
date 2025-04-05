"""
PandaMap-Color: Protein-Ligand Interaction Mapper with customizable color schemes
"""

__version__ = "2.1"
__author__ = "Pritam Kumar Panda"

# Import main components so they're available directly from the package
import warnings
from importlib.metadata import version
import threading
from .pandamap import PandaMapColor
from .colorschemes import COLOR_SCHEMES
from .visualization import visualize
from .ligand import LigandStructure

# For backward compatibility
ProtLigMapper = PandaMapColor

# Attach the visualization function to the PandaMapColor class
PandaMapColor.visualize = visualize

__all__ = [
    'PandaMapColor',
    'ProtLigMapper',  # Compatibility alias
    'LigandStructure',
    'COLOR_SCHEMES'
]
__version__ = "2.1"  # Keep this as fallback if importlib fails

# --- Auto-update checker (non-blocking) ---
def _check_for_updates():
    """Check PyPI for newer versions and notify user with red box."""
    try:
        import requests  # Lazy import
        import importlib.util
        import shutil

        package_name = "pandamap-color"
        current_version = version(package_name)

        response = requests.get(
            f"https://pypi.org/pypi/{package_name}/json",
            timeout=2
        )
        latest_version = response.json()["info"]["version"]

        if current_version != latest_version:
            message = (
                f"\n🚨 [bold red]PandaMap-Color {latest_version} is available![/bold red] "
                f"[dim](you have {current_version})[/dim]\n\n"
                f"[yellow]Update with:[/yellow] [green]pip install --upgrade {package_name}[/green]\n"
                f"[dim]To disable update checks, set: PANDAMAP-COLOR_NO_UPDATE_CHECK=1[/dim]\n"
            )

            # Use rich if available and stdout is a terminal
            if shutil.which("rich"):
                try:
                    from rich.console import Console
                    console = Console()
                    console.print(message)
                except ImportError:
                    warnings.warn(message, UserWarning, stacklevel=2)
            else:
                warnings.warn(message, UserWarning, stacklevel=2)

    except Exception:
        pass  # Don't crash anything if this fails

# Run check only if not disabled
if not __import__('os').getenv("PANDAMAP-COLOR_NO_UPDATE_CHECK"):
    threading.Thread(target=_check_for_updates, daemon=True).start()
