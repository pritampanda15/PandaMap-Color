"""
PandaMap-Color: Protein-Ligand Interaction Mapper with customizable color schemes
"""

__version__ = "2.0.0"
__author__ = "Pritam Kumar Panda"

# Import main components so they're available directly from the package
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
