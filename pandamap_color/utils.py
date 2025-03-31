
"""
Utility functions for PandaMap-Color.
"""

import os
import tempfile

# BioPython imports
from Bio.PDB import PDBParser, NeighborSearch, Selection

# Define three_to_one conversion manually if import isn't available
try:
    from Bio.PDB.Polypeptide import three_to_one
except ImportError:
    # Define the conversion dictionary manually
    _aa_index = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    def three_to_one(residue):
        """Convert amino acid three letter code to one letter code."""
        if residue in _aa_index:
            return _aa_index[residue]
        else:
            return "X"  # Unknown amino acid

    def create_rounded_rectangle(self, xy, width, height, radius=0.1, **kwargs):
        """Create a rectangle with rounded corners using a more compatible approach."""
        # Create vertices for the rounded rectangle
        x, y = xy
    
        # Create a Path with rounded corners
        verts = [
            (x + radius, y),                      # Start
            (x + width - radius, y),              # Top edge
            (x + width, y),                       # Top-right curve start
            (x + width, y + radius),              # Right edge start
            (x + width, y + height - radius),     # Right edge
            (x + width, y + height),              # Bottom-right curve start
            (x + width - radius, y + height),     # Bottom edge start
            (x + radius, y + height),             # Bottom edge
            (x, y + height),                      # Bottom-left curve start
            (x, y + height - radius),             # Left edge start
            (x, y + radius),                      # Left edge
            (x, y),                               # Top-left curve start
            (x + radius, y),                      # Back to start
        ]
    
        # Add to plot as a Polygon instead of trying to use BoxStyle
        rect = Polygon(verts, closed=True, **kwargs)
        return rect