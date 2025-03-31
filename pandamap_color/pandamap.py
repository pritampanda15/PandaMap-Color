
"""
Main PandaMap-Color class for analyzing protein-ligand interactions.
"""

import os
from collections import defaultdict
from Bio.PDB import PDBParser, NeighborSearch
import sys
import math
import argparse
import json
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

from .colorschemes import COLOR_SCHEMES, load_custom_color_scheme
from .ligand import LigandStructure

class PandaMapColor:
    """
    Class for analyzing protein-ligand interactions and creating 
    visualizations with customizable styling.
    """
    
    def __init__(self, pdb_file, ligand_resname=None, color_scheme='default', use_enhanced_styling=True):
        """
        Initialize with a PDB file containing a protein-ligand complex.
        
        Parameters:
        -----------
        pdb_file : str
            Path to the PDB file with protein and ligand
        ligand_resname : str, optional
            Specific residue name of the ligand to focus on
        color_scheme : str or dict
            Color scheme to use, either a key from COLOR_SCHEMES or a custom dictionary
        use_enhanced_styling : bool
            Whether to use enhanced styling effects (gradients, shadows, etc.)
        """
        self.pdb_file = pdb_file
        self.ligand_resname = ligand_resname
        self.use_enhanced_styling = use_enhanced_styling
        
        # Set color scheme
        if isinstance(color_scheme, str):
            if color_scheme in COLOR_SCHEMES:
                self.colors = COLOR_SCHEMES[color_scheme]
            else:
                print(f"Warning: Unknown color scheme '{color_scheme}'. Using default.")
                self.colors = COLOR_SCHEMES['default']
        elif isinstance(color_scheme, dict):
            # Merge with default scheme to ensure all required colors are present
            self.colors = COLOR_SCHEMES['default'].copy()
            self.colors.update(color_scheme)
        else:
            self.colors = COLOR_SCHEMES['default']
        
        # Parse the PDB file
        try:
            self.parser = PDBParser(QUIET=True)
            self.structure = self.parser.get_structure('complex', pdb_file)
            self.model = self.structure[0]
        except Exception as e:
            raise ValueError(f"Failed to parse PDB file: {e}")
        
        # Separate ligand from protein
        self.protein_atoms = []
        self.ligand_atoms = []
        self.protein_residues = {}
        self.ligand_residue = None
        
        for residue in self.model.get_residues():
            # Store ligand atoms (HETATM records)
            if residue.id[0] != ' ':  # Non-standard residue (HETATM)
                if ligand_resname is None or residue.resname == ligand_resname:
                    for atom in residue:
                        self.ligand_atoms.append(atom)
                    if self.ligand_residue is None:
                        self.ligand_residue = residue
            else:  # Standard residues (protein)
                res_id = (residue.resname, residue.id[1])
                self.protein_residues[res_id] = residue
                for atom in residue:
                    self.protein_atoms.append(atom)
        
        # Check if we found a ligand
        if not self.ligand_atoms:
            raise ValueError("No ligand (HETATM) found in the PDB file.")
        
        # Storage for the interaction data
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': []
        }
        
        # Will store residues that interact with the ligand
        self.interacting_residues = set()
        
        # For solvent accessibility information (simplified)
        self.solvent_accessible = set()
        
        # Create the ligand structure
        self.ligand_structure = LigandStructure(
            self.ligand_atoms, 
            color_scheme=color_scheme, 
            use_enhanced_styling=use_enhanced_styling
        )
        
    def detect_interactions(self, 
                           h_bond_cutoff=3.5, 
                           pi_stack_cutoff=5.5,
                           hydrophobic_cutoff=4.0):
        """
        Detect all interactions between protein and ligand.
        
        Parameters:
        -----------
        h_bond_cutoff : float
            Distance cutoff for hydrogen bonds in Angstroms
        pi_stack_cutoff : float
            Distance cutoff for pi-stacking interactions in Angstroms
        hydrophobic_cutoff : float
            Distance cutoff for hydrophobic interactions in Angstroms
        """
        # Use neighbor search for efficiency
        ns = NeighborSearch(self.protein_atoms)
        max_cutoff = max(h_bond_cutoff, pi_stack_cutoff, hydrophobic_cutoff)
        
        # Define amino acid categories
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        h_bond_donors = {'ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'TRP'}
        h_bond_acceptors = {'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        neg_charged = {'ASP', 'GLU'}
        amide_residues = {'ASN', 'GLN'}
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        
        # Store all potential interactions first
        all_interactions = {
            interaction_type: []
            for interaction_type in self.interactions.keys()
        }
        
        # Check each ligand atom for interactions
        for lig_atom in self.ligand_atoms:
            # Find protein atoms within cutoff distance
            nearby_atoms = ns.search(lig_atom.get_coord(), max_cutoff)
            
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = lig_atom - prot_atom
                
                # Store interacting residue for later visualization
                res_id = (prot_res.resname, prot_res.id[1])
                self.interacting_residues.add(res_id)
                
                # Create a common structure for the interaction info
                interaction_info = {
                    'ligand_atom': lig_atom,
                    'protein_atom': prot_atom,
                    'protein_residue': prot_res,
                    'distance': distance,
                    'res_id': res_id
                }
                
                # Check hydrogen bonds
                if (distance <= h_bond_cutoff and 
                    lig_atom.element in ['N', 'O'] and prot_atom.element in ['N', 'O']):
                    all_interactions['hydrogen_bonds'].append(interaction_info)
                
                # Check pi-pi stacking
                if (distance <= pi_stack_cutoff and 
                    prot_res.resname in aromatic_residues and 
                    lig_atom.element == 'C' and prot_atom.element == 'C'):
                    all_interactions['pi_pi_stacking'].append(interaction_info)
                
                # Check carbon-pi interactions
                if (distance <= pi_stack_cutoff and 
                    prot_res.resname in aromatic_residues and 
                    lig_atom.element == 'C'):
                    all_interactions['carbon_pi'].append(interaction_info)
                
                # Check donor-pi interactions
                if (distance <= pi_stack_cutoff and 
                    prot_res.resname in neg_charged and 
                    lig_atom.element == 'C'):
                    all_interactions['donor_pi'].append(interaction_info)
                
                # Check amide-pi interactions
                if (distance <= pi_stack_cutoff and 
                    prot_res.resname in amide_residues and 
                    lig_atom.element == 'C'):
                    all_interactions['amide_pi'].append(interaction_info)
                
                # Check hydrophobic interactions
                if (distance <= hydrophobic_cutoff and 
                    prot_res.resname in hydrophobic_residues and 
                    lig_atom.element == 'C' and prot_atom.element == 'C'):
                    all_interactions['hydrophobic'].append(interaction_info)
        
        # Deduplicate interactions - keep only the closest interaction of each type per residue
        for interaction_type, interactions_list in all_interactions.items():
            # Group interactions by residue
            by_residue = defaultdict(list)
            for interaction in interactions_list:
                res_id = interaction['res_id']
                by_residue[res_id].append(interaction)
            
            # Keep only the closest interaction per residue
            self.interactions[interaction_type] = []
            for res_id, res_interactions in by_residue.items():
                if res_interactions:
                    closest = min(res_interactions, key=lambda x: x['distance'])
                    self.interactions[interaction_type].append(closest)
    
    def estimate_solvent_accessibility(self):
        """
        Estimate which residues might be solvent accessible.
        This is a simplified approach since we're trying to match the example image.
        In a real implementation, you'd use DSSP or a similar tool.
        """
        # For simplicity, mark all residues as solvent accessible
        # In a real implementation, you'd use a proper algorithm
        self.solvent_accessible = self.interacting_residues.copy()
    
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
    
        
    def run_analysis(self, output_file=None):
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved. If None, a default name will be generated.
            
        Returns:
        --------
        str : Path to the generated visualization file or error message
        """
        try:
            # Ensure we have a valid output filename
            if output_file is None:
                base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
                output_file = f"{base_name}_interactions.png"
                print(f"Using default output filename: {output_file}")
            
            # Detect protein-ligand interactions
            print("Detecting interactions...")
            self.detect_interactions()
            
            # Estimate solvent accessibility
            print("Estimating solvent accessibility...")
            self.estimate_solvent_accessibility()
            
            # Generate visualization
            print("Generating visualization...")
            result = self.visualize(output_file=output_file)
            
            return result
            
        except Exception as e:
            error_msg = f"Error in analysis pipeline: {e}"
            print(f"✗ {error_msg}")
            return error_msg