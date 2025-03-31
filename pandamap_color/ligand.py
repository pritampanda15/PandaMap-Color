"""
Ligand structure handling for PandaMap-Color.

This module contains the LigandStructure class for creating
2D representations of ligand molecules.
"""
import os
import sys
import math
import argparse
import json
from collections import defaultdict
import matplotlib as mpl

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Rectangle, Polygon, PathPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.font_manager as fm
from matplotlib.path import Path
from matplotlib.collections import LineCollection

# Import COLOR_SCHEMES from colorschemes.py
from .colorschemes import COLOR_SCHEMES

class LigandStructure:
    """
    Class to create a 2D representation of a ligand structure
    with customizable aesthetics.
    """
    
    def __init__(self, ligand_atoms, color_scheme='default', use_enhanced_styling=True):
        """
        Initialize with a list of ligand atoms from a BioPython structure.
        
        Parameters:
        -----------
        ligand_atoms : list
            List of BioPython Atom objects from the ligand
        color_scheme : str or dict
            Color scheme to use, either a key from COLOR_SCHEMES or a custom dictionary
        use_enhanced_styling : bool
            Whether to use enhanced styling effects (gradients, shadows, etc.)
        """
        self.ligand_atoms = ligand_atoms
        self.atom_coords = {}
        self.use_enhanced_styling = use_enhanced_styling
        
        # Set color scheme
        if isinstance(color_scheme, str):
            if color_scheme in COLOR_SCHEMES:
                self.element_colors = COLOR_SCHEMES[color_scheme]['element_colors']
            else:
                print(f"Warning: Unknown color scheme '{color_scheme}'. Using default.")
                self.element_colors = COLOR_SCHEMES['default']['element_colors']
        elif isinstance(color_scheme, dict) and 'element_colors' in color_scheme:
            self.element_colors = color_scheme['element_colors']
        else:
            self.element_colors = COLOR_SCHEMES['default']['element_colors']
        
        # Record atom coordinates and elements
        for atom in ligand_atoms:
            atom_id = atom.get_id()
            self.atom_coords[atom_id] = {
                'element': atom.element,
                'coord': atom.get_coord(),  # 3D coordinates from PDB
                'name': atom.get_name()     # Atom name for labeling
            }
    
    def generate_2d_coords(self):
        """
        Generate 2D coordinates for the ligand atoms based on their 3D coordinates.
        Uses PCA to find the best 2D projection plane.
        
        Returns:
        --------
        dict : Dictionary mapping atom IDs to 2D coordinates
        """
        if not self.atom_coords:
            return {}
            
        # Get all 3D coordinates and find center
        all_coords = np.array([info['coord'] for info in self.atom_coords.values()])
        center = np.mean(all_coords, axis=0)
        
        # Subtract center to center the molecule
        centered_coords = all_coords - center
        
        # PCA-like approach to find main plane
        cov_matrix = np.cov(centered_coords.T)
        
        try:
            # Get eigenvectors and eigenvalues
            eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
            
            # Sort by eigenvalue in descending order
            idx = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            
            # Use the first two eigenvectors to define the plane
            plane_vectors = eigenvectors[:, :2]
            
            # Project the centered coordinates onto the plane
            projected_coords = np.dot(centered_coords, plane_vectors)
            
            # Scale to fit nicely in the visualization
            max_dim = np.max(np.abs(projected_coords))
            scaling_factor = 50.0 / max_dim if max_dim > 0 else 1.0
            projected_coords *= scaling_factor
            
            # Store the 2D coordinates
            coords_2d = {}
            for i, atom_id in enumerate(self.atom_coords.keys()):
                coords_2d[atom_id] = projected_coords[i]
                
        except np.linalg.LinAlgError:
            # Fallback if eigendecomposition fails
            print("Warning: Could not compute optimal projection. Using simple XY projection.")
            coords_2d = {}
            for atom_id, info in self.atom_coords.items():
                # Simple scaling of x, y coordinates
                coords_2d[atom_id] = np.array([info['coord'][0], info['coord'][1]]) * 10.0
        
        return coords_2d
    
    def find_bonds(self, distance_threshold=2.0):
        """
        Find bonds between atoms based on distance with improved bond detection.
        
        Parameters:
        -----------
        distance_threshold : float
            Maximum distance between atoms to be considered bonded (in Angstroms)
            
        Returns:
        --------
        list : List of tuples (atom_id1, atom_id2, bond_type) representing bonds
        """
        bonds = []
        atom_ids = list(self.atom_coords.keys())
        
        # Define common bond distances for better accuracy
        typical_bond_lengths = {
            ('C', 'C'): 1.5,    # Single bond
            ('C', 'N'): 1.4,    # Single bond
            ('C', 'O'): 1.4,    # Single bond
            ('C', 'S'): 1.8,    # Single bond
            ('N', 'N'): 1.4,    # Single bond
            ('N', 'O'): 1.4,    # Single bond
            ('O', 'P'): 1.6,    # Single bond
            ('S', 'O'): 1.5,    # Single bond
        }
        
        # Tolerance factor for bond detection
        tolerance = 0.45  # Angstroms
        
        for i in range(len(atom_ids)):
            for j in range(i+1, len(atom_ids)):
                atom1_id = atom_ids[i]
                atom2_id = atom_ids[j]
                
                elem1 = self.atom_coords[atom1_id]['element']
                elem2 = self.atom_coords[atom2_id]['element']
                
                coord1 = self.atom_coords[atom1_id]['coord']
                coord2 = self.atom_coords[atom2_id]['coord']
                
                # Calculate Euclidean distance
                distance = np.linalg.norm(coord1 - coord2)
                
                # Get typical bond length for this pair of elements
                elems = tuple(sorted([elem1, elem2]))
                typical_length = typical_bond_lengths.get(elems, 1.5)  # Default if not in dictionary
                
                # Determine if atoms are bonded
                if distance < (typical_length + tolerance):
                    # Simple estimation of bond type based on distance
                    bond_type = 'single'  # Default is single bond
                    
                    # For now, we'll just use single bonds for simplicity
                    # In a more sophisticated implementation, you'd determine
                    # double and triple bonds based on geometric and electronic considerations
                    
                    bonds.append((atom1_id, atom2_id, bond_type))
        
        return bonds
    
    def draw_on_axes(self, ax, center=(0, 0), radius=80):
        """
        Draw a 2D representation of the ligand with customizable styling.
        
        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw
        center : tuple
            The (x, y) coordinates where the center of the molecule should be
        radius : float
            The approximate radius the molecule should occupy
            
        Returns:
        --------
        dict : Dictionary mapping atom IDs to their 2D positions in the plot
        """
        # Generate 2D coordinates
        coords_2d = self.generate_2d_coords()
        
        if not coords_2d:
            # If we couldn't generate coordinates, draw a simple placeholder
            print("Warning: Could not generate ligand coordinates. Drawing placeholder.")
            circle = Circle(center, radius/2, fill=False, edgecolor='black', linestyle='-')
            ax.add_patch(circle)
            ax.text(center[0], center[1], "Ligand", ha='center', va='center')
            return {}
            
        # Find bonds
        bonds = self.find_bonds()
        
        # Scale coordinates to fit within the specified radius
        all_coords = np.array(list(coords_2d.values()))
        max_extent = np.max(np.abs(all_coords))
        scaling_factor = radius / (max_extent * 1.2)  # Leave some margin
        
        # Create a mapping of atom IDs to positions in the plot
        atom_positions = {}
        
        # Draw bonds first (so they're below atoms)
        for atom1_id, atom2_id, bond_type in bonds:
            pos1 = coords_2d[atom1_id] * scaling_factor + center
            pos2 = coords_2d[atom2_id] * scaling_factor + center
            
            # Get elements for determining bond appearance
            elem1 = self.atom_coords[atom1_id]['element']
            elem2 = self.atom_coords[atom2_id]['element']
            
            if self.use_enhanced_styling:
                # Draw bond with gradient color from one atom to the other
                color1 = self.element_colors.get(elem1, 'gray')
                color2 = self.element_colors.get(elem2, 'gray')
                
                # Create a custom colormap for this bond
                cmap = LinearSegmentedColormap.from_list("bond_cmap", [color1, color2])
                
                # Draw single bond as a gradient line
                x = np.linspace(pos1[0], pos2[0], 100)
                y = np.linspace(pos1[1], pos2[1], 100)
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                
                # Create a line collection with the gradient colormap
                lc = LineCollection(segments, cmap=cmap, linewidth=2.5, alpha=0.8, zorder=2)
                
                # Set color values for the gradient
                lc.set_array(np.linspace(0, 1, len(segments)))
                ax.add_collection(lc)
            else:
                # Simple bond drawing without gradient
                line = Line2D([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                             color='black', linewidth=1.5, zorder=2)
                ax.add_line(line)
        
        # Draw atoms as circles with styling based on settings
        for atom_id, coord in coords_2d.items():
            # Scale and shift the position
            pos = coord * scaling_factor + center
            atom_positions[atom_id] = pos
            
            element = self.atom_coords[atom_id]['element']
            color = self.element_colors.get(element, 'gray')
            
            # Determine size based on element (larger for heavier atoms)
            size = 9 if element in ['C'] else 12
            
            if self.use_enhanced_styling:
                # Draw atom with a slight shadow for 3D effect
                shadow = Circle(pos + np.array([1, 1]), size, 
                               facecolor='black', alpha=0.2, zorder=2.8)
                ax.add_patch(shadow)
                
                # Main atom circle with enhanced styling
                circle = Circle(pos, size, facecolor=color, edgecolor='black', 
                               linewidth=1, alpha=0.9, zorder=3)
            else:
                # Simpler atom drawing without shadow
                circle = Circle(pos, size, facecolor=color, edgecolor='black', 
                               linewidth=1, alpha=0.8, zorder=3)
            
            ax.add_patch(circle)
            
            # Add element label (except for carbon)
            if element != 'C':
                text = ax.text(pos[0], pos[1], element, ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='white', zorder=4)
                
                # Add a subtle shadow to the text for better visibility if using enhanced styling
                if self.use_enhanced_styling:
                    text.set_path_effects([path_effects.withStroke(linewidth=1.5, 
                                                                foreground='black')])
        
        return atom_positions