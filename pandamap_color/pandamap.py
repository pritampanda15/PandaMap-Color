"""
Main PandaMap-Color class for analyzing protein-ligand interactions.
Enhanced with improved interaction detection from core.py and improved_interaction_detection.py
"""

import os
import sys
import math
import numpy as np
import tempfile
from collections import defaultdict, namedtuple
from Bio.PDB import PDBParser, NeighborSearch, Selection, PDBIO
import matplotlib as mpl
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

# Constants for interaction detection thresholds based on PLIP
DISTANCE_CUTOFFS = {
    'hydrogen_bond': 3.5,       # Max distance for hydrogen bonds
    'hydrophobic': 4.0,         # Max distance for hydrophobic interactions
    'pi_stacking': 5.5,         # Max distance for pi-stacking
    'pi_cation': 6.0,           # Max distance for pi-cation interactions
    'salt_bridge': 5.5,         # Max distance for salt bridges
    'halogen': 3.5,             # Max distance for halogen bonds
    'water_bridge': 3.5,        # Max distance for water bridges
    'min_dist': 1.5             # Minimum distance for any interaction
}

ANGLE_CUTOFFS = {
    'pi_stack_parallel': 30.0,      # Max deviation from parallel planes (degrees)
    'pi_stack_perpendicular': 30.0, # Max deviation from perpendicular planes (degrees)
    'pi_cation': 30.0,              # Max deviation for pi-cation (degrees)
    'hbond_donor': 120.0,           # Min angle for hydrogen bond donor (degrees)
    'halogen_donor': 140.0,         # Optimal angle for halogen donor (C-X...Y)
    'halogen_acceptor': 120.0,      # Optimal angle for halogen acceptor (X...Y-C)
    'halogen_angle_dev': 30.0       # Max deviation from optimal halogen bond angles
}

class PandaMapColor:
    """
    Class for analyzing protein-ligand interactions and creating 
    visualizations with customizable styling.
    Enhanced with improved interaction detection algorithms.
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
        
        # Track metal and halogen atoms for extended interaction detection
        self.metal_atoms = []
        self.halogen_atoms = []
        
        for residue in self.model.get_residues():
            # Store ligand atoms (HETATM records)
            if residue.id[0] != ' ':  # Non-standard residue (HETATM)
                if ligand_resname is None or residue.resname == ligand_resname:
                    for atom in residue:
                        self.ligand_atoms.append(atom)
                        # Identify halogen atoms in ligand
                        if atom.element in ['F', 'Cl', 'Br', 'I']:
                            self.halogen_atoms.append(atom)
                    if self.ligand_residue is None:
                        self.ligand_residue = residue
                else:
                    # Check for metal ions (typically single-atom HETATM residues)
                    if len(list(residue.get_atoms())) == 1:
                        atom = next(residue.get_atoms())
                        if atom.element in ['MG', 'ZN', 'CA', 'FE', 'MN', 'CU', 'NA', 'K', 'LI', 'CO', 'NI']:
                            self.metal_atoms.append(atom)
            else:  # Standard residues (protein)
                res_id = (residue.resname, residue.id[1])
                self.protein_residues[res_id] = residue
                for atom in residue:
                    self.protein_atoms.append(atom)
        
        # Check if we found a ligand
        if not self.ligand_atoms:
            raise ValueError("No ligand (HETATM) found in the PDB file.")
        
        # Storage for the interaction data - EXPANDED with additional interaction types
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': [],
            'ionic': [],           # Ionic interactions
            'halogen_bonds': [],   # Halogen bonds
            'cation_pi': [],       # Cation-pi interactions
            'metal_coordination': [], # Metal coordination
            'salt_bridge': [],     # Salt bridges
            'covalent': [],        # Covalent bonds
            'alkyl_pi': [],        # Alkyl-Pi interactions
            'attractive_charge': [], # Attractive charge interactions
            'pi_cation': [],       # Pi-cation interactions
            'repulsion': []        # Repulsion interactions
        }
        
        # Will store residues that interact with the ligand
        self.interacting_residues = set()
        
        # Store interaction directionality (protein->ligand, ligand->protein, or both)
        self.interaction_direction = {}
        
        # For solvent accessibility information
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
                           hydrophobic_cutoff=4.0,
                           ionic_cutoff=4.0,
                           halogen_bond_cutoff=3.5,
                           metal_coord_cutoff=2.8,
                           covalent_cutoff=2.1):
        """
        Detect all interactions between protein and ligand with enhanced detection.
        
        Parameters:
        -----------
        h_bond_cutoff : float
            Distance cutoff for hydrogen bonds in Angstroms
        pi_stack_cutoff : float
            Distance cutoff for pi-stacking interactions in Angstroms
        hydrophobic_cutoff : float
            Distance cutoff for hydrophobic interactions in Angstroms
        ionic_cutoff : float
            Distance cutoff for ionic/salt bridge interactions in Angstroms
        halogen_bond_cutoff : float
            Distance cutoff for halogen bonds in Angstroms
        metal_coord_cutoff : float
            Distance cutoff for metal coordination in Angstroms
        covalent_cutoff : float
            Distance cutoff for covalent bonds in Angstroms
        """
        print("Detecting interactions with enhanced criteria...")
        # Use neighbor search for efficiency
        ns = NeighborSearch(self.protein_atoms)
        max_cutoff = max(h_bond_cutoff, pi_stack_cutoff, hydrophobic_cutoff, 
                      ionic_cutoff, halogen_bond_cutoff, metal_coord_cutoff,
                      covalent_cutoff)
        
        # Define amino acid categories
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        h_bond_donors = {'ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'TRP'}
        h_bond_acceptors = {'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        neg_charged = {'ASP', 'GLU'}
        pos_charged = {'ARG', 'LYS', 'HIS'}
        amide_residues = {'ASN', 'GLN'}
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        alkyl_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'}
        
        # Reset the interaction data
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': [],
            'ionic': [],
            'halogen_bonds': [],
            'cation_pi': [],
            'metal_coordination': [],
            'salt_bridge': [],
            'covalent': [],
            'alkyl_pi': [],
            'attractive_charge': [],
            'pi_cation': [],
            'repulsion': []
        }
        
        # Reset interacting residues
        self.interacting_residues = set()
        
        # Helper functions for determining atom types
        def is_halogen_acceptor(atom):
            return atom.element in ['O', 'N', 'S'] or (atom.element == 'C' and atom.name in ['CE1', 'CD2'])
        
        def is_metal_coordinator(atom):
            return atom.element in ['O', 'N', 'S'] or atom.name in ['SD', 'OD1', 'OD2', 'OE1', 'OE2', 'NE2', 'ND1']
        
        # Initialize temporary storage for all interactions
        all_interactions = {key: [] for key in self.interactions.keys()}
        
        # Check each ligand atom for interactions
        for lig_atom in self.ligand_atoms:
            # Find protein atoms within cutoff distance
            nearby_atoms = ns.search(lig_atom.get_coord(), max_cutoff)
            
            # Track ionizable ligand atoms for charged interactions
            is_lig_pos_charged = lig_atom.element == 'N' and not any(
                a.element == 'C' for a in self.ligand_atoms 
                if a.element != 'H' and np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 1.6
            )
            is_lig_neg_charged = lig_atom.element == 'O' and not any(
                a.element == 'C' for a in self.ligand_atoms 
                if a.element != 'H' and np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 1.6
            )
            
            # Track aromatic ligand characteristics
            # An atom is potentially part of an aromatic ring if it's a carbon 
            # with at least 2 other carbon neighbors at aromatic bond distances
            is_lig_aromatic = lig_atom.element == 'C' and len([
                a for a in self.ligand_atoms 
                if a.element == 'C' and 1.2 < np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 2.8
            ]) >= 2
            
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = lig_atom - prot_atom
                
                # Skip if distance is too small (likely a clash or error)
                if distance < 1.5:
                    continue
                
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
                
                # 1. Hydrogen bonds - N and O atoms within cutoff
                if (distance <= h_bond_cutoff and distance >= 2.4 and  # Added minimum distance
                    lig_atom.element in ['N', 'O'] and prot_atom.element in ['N', 'O']):
                    
                    all_interactions['hydrogen_bonds'].append(interaction_info)
                    
                    # Determine directionality
                    interaction_key = (res_id, 'hydrogen_bonds')
                    is_donor_prot = prot_res.resname in h_bond_donors and prot_atom.element == 'N'
                    is_acceptor_prot = prot_res.resname in h_bond_acceptors and prot_atom.element in ['O', 'N']
                    is_donor_lig = lig_atom.element == 'N'
                    is_acceptor_lig = lig_atom.element in ['O', 'N']
                    
                    if (is_donor_prot and is_acceptor_lig) and (is_donor_lig and is_acceptor_prot):
                        self.interaction_direction[interaction_key] = 'bidirectional'
                    elif is_donor_prot and is_acceptor_lig:
                        self.interaction_direction[interaction_key] = 'protein_to_ligand'
                    elif is_donor_lig and is_acceptor_prot:
                        self.interaction_direction[interaction_key] = 'ligand_to_protein'
                    else:
                        self.interaction_direction[interaction_key] = 'bidirectional'
                
                # 2. Pi-cation interactions
                if distance <= pi_stack_cutoff:
                    is_prot_pos_charged = prot_res.resname in pos_charged and prot_atom.element in ['N']
                    is_lig_aromatic_ring = is_lig_aromatic and lig_atom.element == 'C'
                    
                    if is_prot_pos_charged and is_lig_aromatic_ring:
                        all_interactions['pi_cation'].append(interaction_info)
                
                # 3. Alkyl-Pi interactions
                if distance <= pi_stack_cutoff:
                    # Check for alkyl groups in protein interacting with aromatic ligand
                    is_prot_alkyl = prot_res.resname in alkyl_residues and prot_atom.element == 'C'
                    
                    if is_prot_alkyl and is_lig_aromatic:
                        all_interactions['alkyl_pi'].append(interaction_info)
                    
                    # Check for aromatic groups in protein interacting with alkyl in ligand
                    is_prot_aromatic = prot_res.resname in aromatic_residues and prot_atom.element == 'C'
                    is_lig_alkyl = lig_atom.element == 'C' and not is_lig_aromatic
                    
                    if is_prot_aromatic and is_lig_alkyl:
                        all_interactions['alkyl_pi'].append(interaction_info)
                
                # 4. Attractive charge interactions
                if distance <= ionic_cutoff:
                    is_prot_pos = prot_res.resname in pos_charged and prot_atom.element == 'N'
                    is_prot_neg = prot_res.resname in neg_charged and prot_atom.element == 'O'
                    
                    if (is_prot_pos and is_lig_neg_charged) or (is_prot_neg and is_lig_pos_charged):
                        all_interactions['attractive_charge'].append(interaction_info)
                
                # 5. Repulsion interactions
                if distance <= ionic_cutoff * 1.5:  # Larger cutoff for repulsion
                    is_prot_pos_charged = prot_res.resname in pos_charged and prot_atom.element in ['N']
                    is_prot_neg_charged = prot_res.resname in neg_charged and prot_atom.element in ['O']
                    
                    if (is_prot_pos_charged and is_lig_pos_charged) or (is_prot_neg_charged and is_lig_neg_charged):
                        all_interactions['repulsion'].append(interaction_info)
                
                # 6. Pi-pi stacking
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in aromatic_residues and is_lig_aromatic and lig_atom.element == 'C' and prot_atom.element == 'C':
                        all_interactions['pi_pi_stacking'].append(interaction_info)
                
                # 7. Carbon-Pi interactions
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in aromatic_residues and lig_atom.element == 'C':
                        all_interactions['carbon_pi'].append(interaction_info)
                
                # 8. Donor-Pi interactions
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in neg_charged and is_lig_aromatic:
                        all_interactions['donor_pi'].append(interaction_info)
                
                # 9. Amide-Pi interactions
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in amide_residues and is_lig_aromatic:
                        all_interactions['amide_pi'].append(interaction_info)
                
                # 10. Hydrophobic interactions
                if distance <= hydrophobic_cutoff:
                    if (prot_res.resname in hydrophobic_residues and 
                        lig_atom.element == 'C' and prot_atom.element == 'C'):
                        all_interactions['hydrophobic'].append(interaction_info)
                
                # 11. Ionic/salt bridge interactions
                if distance <= ionic_cutoff:
                    is_prot_pos = prot_res.resname in pos_charged and prot_atom.element == 'N'
                    is_prot_neg = prot_res.resname in neg_charged and prot_atom.element == 'O'
                    
                    if (is_prot_pos and is_lig_neg_charged) or (is_prot_neg and is_lig_pos_charged):
                        all_interactions['ionic'].append(interaction_info)
                        all_interactions['salt_bridge'].append(interaction_info)
                
                # 12. Cation-Pi interactions
                if distance <= pi_stack_cutoff:
                    if ((prot_res.resname in pos_charged and is_lig_aromatic) or 
                        (prot_res.resname in aromatic_residues and is_lig_pos_charged)):
                        all_interactions['cation_pi'].append(interaction_info)
                
                # 13. Covalent bonds
                if distance <= covalent_cutoff:
                    # Only include likely covalent bonds involving specific residues and atom types
                    if ((prot_res.resname == 'CYS' and prot_atom.name == 'SG') or 
                        (prot_res.resname == 'SER' and prot_atom.name == 'OG') or
                        (prot_res.resname == 'LYS' and prot_atom.name == 'NZ') or
                        (prot_res.resname == 'HIS' and prot_atom.name in ['ND1', 'NE2'])):
                        all_interactions['covalent'].append(interaction_info)
        
        # Handle halogen bonds separately
        for halogen_atom in self.halogen_atoms:
            nearby_atoms = ns.search(halogen_atom.get_coord(), halogen_bond_cutoff)
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = halogen_atom - prot_atom
                
                if 1.5 < distance <= halogen_bond_cutoff and is_halogen_acceptor(prot_atom):
                    res_id = (prot_res.resname, prot_res.id[1])
                    self.interacting_residues.add(res_id)
                    
                    interaction_info = {
                        'ligand_atom': halogen_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance,
                        'res_id': res_id
                    }
                    all_interactions['halogen_bonds'].append(interaction_info)
        
        # Handle metal coordination
        for metal_atom in self.metal_atoms:
            nearby_atoms = ns.search(metal_atom.get_coord(), metal_coord_cutoff)
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = metal_atom - prot_atom
                
                if distance <= metal_coord_cutoff and is_metal_coordinator(prot_atom):
                    res_id = (prot_res.resname, prot_res.id[1])
                    self.interacting_residues.add(res_id)
                    
                    interaction_info = {
                        'ligand_atom': metal_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance,
                        'res_id': res_id
                    }
                    all_interactions['metal_coordination'].append(interaction_info)
        
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
        
        # Apply post-processing filters for better chemical plausibility
        self._filter_interactions_chemically()
        
        # Count total interactions
        total_interactions = sum(len(ints) for ints in self.interactions.values())
        print(f"Detected {total_interactions} interactions across {len(self.interacting_residues)} residues")
        
        return self.interactions

    def _filter_interactions_chemically(self):
        """Apply chemical knowledge to filter out implausible interactions."""
        # Define residue categories for filtering
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        charged_residues = {'ASP', 'GLU', 'LYS', 'ARG', 'HIS'}
        neg_charged = {'ASP', 'GLU'}
        pos_charged = {'LYS', 'ARG', 'HIS'}
        amide_residues = {'ASN', 'GLN'}
        
        # Filter pi-stacking - require at least one aromatic residue
        self.interactions['pi_pi_stacking'] = [
            i for i in self.interactions['pi_pi_stacking']
            if i['protein_residue'].resname in aromatic_residues and i['distance'] < 5.5
        ]
        
        # Filter ionic/salt bridge - require charged residues
        self.interactions['ionic'] = [
            i for i in self.interactions['ionic']
            if i['protein_residue'].resname in charged_residues and i['distance'] < 4.0
        ]
        self.interactions['salt_bridge'] = [
            i for i in self.interactions['salt_bridge']
            if i['protein_residue'].resname in charged_residues and i['distance'] < 4.0
        ]
        
        # Filter covalent bonds - require very close distance
        self.interactions['covalent'] = [
            i for i in self.interactions['covalent']
            if i['distance'] < 2.1
        ]
        
        # Filter amide-pi - require amide residues
        self.interactions['amide_pi'] = [
            i for i in self.interactions['amide_pi']
            if i['protein_residue'].resname in amide_residues and i['distance'] < 5.5
        ]
        
        # Filter donor-pi - require negatively charged residues
        self.interactions['donor_pi'] = [
            i for i in self.interactions['donor_pi']
            if i['protein_residue'].resname in neg_charged and i['distance'] < 5.5
        ]
        
        # Filter carbon-pi - require aromatic residues
        self.interactions['carbon_pi'] = [
            i for i in self.interactions['carbon_pi']
            if i['protein_residue'].resname in aromatic_residues and i['distance'] < 5.5
        ]
        
        # Filter cation-pi - require charged residues or aromatic residues
        self.interactions['cation_pi'] = [
            i for i in self.interactions['cation_pi']
            if (i['protein_residue'].resname in pos_charged or 
                i['protein_residue'].resname in aromatic_residues) and i['distance'] < 6.0
        ]

    def calculate_realistic_solvent_accessibility(self, probe_radius=1.4, exposure_threshold=0.25, max_percent=0.5):
        """
        Realistic solvent accessibility calculation with proper constraints on number of accessible residues.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        exposure_threshold : float
            Threshold ratio for considering a residue solvent accessible (default: 0.25)
        max_percent : float
            Maximum percentage of interacting residues that can be solvent accessible (default: 0.5)
        """
        print("Calculating realistic solvent accessibility...")
        self.solvent_accessible = set()
        
        # Define which residues are typically surface-exposed
        likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        likely_buried = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'CYS', 'PRO'}
        
        # First get all protein atoms (including non-interacting ones)
        all_protein_atoms = []
        for residue in self.model.get_residues():
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    all_protein_atoms.append(atom)
        
        print(f"Total protein atoms: {len(all_protein_atoms)}")
        print(f"Interacting residues to check: {len(self.interacting_residues)}")
        
        # Calculate protein center
        protein_center = np.zeros(3)
        for atom in all_protein_atoms:
            protein_center += atom.get_coord()
        protein_center /= len(all_protein_atoms) if all_protein_atoms else 1
        
        # Store exposure scores for all residues to sort later
        exposure_scores = {}
        
        # For each interacting residue, estimate accessibility
        for res_id in self.interacting_residues:
            residue = self.protein_residues.get(res_id)
            if residue is None:
                continue
            
            # Calculate residue center
            residue_center = np.zeros(3)
            residue_atoms = list(residue.get_atoms())
            for atom in residue_atoms:
                residue_center += atom.get_coord()
            residue_center /= len(residue_atoms) if residue_atoms else 1
            
            # Calculate vector from protein center to residue center
            direction = residue_center - protein_center
            direction_length = np.linalg.norm(direction)
            if direction_length > 0:
                direction = direction / direction_length
            
            # Count exposed atoms
            exposed_atoms = 0
            total_atoms = len(residue_atoms)
            
            for atom in residue_atoms:
                atom_coord = atom.get_coord()
                
                # Check if atom is on the protein surface
                is_exposed = True
                nearby_atom_count = 0
                
                for other_atom in all_protein_atoms:
                    if other_atom.get_parent() == residue:
                        continue  # Skip atoms in same residue
                    
                    distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                    
                    # Simple distance threshold for exposure
                    if distance < 3.0:
                        nearby_atom_count += 1
                    
                    # If more than 8 atoms are nearby (instead of 12), consider it buried
                    if nearby_atom_count > 8:
                        is_exposed = False
                        break
                
                if is_exposed:
                    exposed_atoms += 1
            
            # Calculate exposure ratio
            exposure_ratio = exposed_atoms / total_atoms if total_atoms > 0 else 0
            
            # Calculate more strict surface bias
            if residue.resname in likely_exposed:
                surface_bias = 1.5  # More bias for typically exposed residues
            elif residue.resname in likely_buried:
                surface_bias = 0.6  # Much lower bias for typically buried residues
            else:
                surface_bias = 1.0
                
            # Calculate distance from surface bias
            # Residues further from center are more likely exposed
            distance_from_center_ratio = min(direction_length / 15.0, 1.0)
            
            # Combined score for exposure - more weight on actual exposure ratio
            exposure_score = (exposure_ratio * 2.0 + surface_bias + distance_from_center_ratio) / 4.0
            
            # Store score for later ranking
            exposure_scores[res_id] = exposure_score
            
            # Only add highest scoring residues directly
            if exposure_score > exposure_threshold:
                self.solvent_accessible.add(res_id)
        
        # Calculate constraints on number of solvent accessible residues
        min_expected = max(1, int(len(self.interacting_residues) * 0.1))  # At least 10%
        max_expected = min(int(len(self.interacting_residues) * max_percent), 
                         len(self.interacting_residues) - 1)  # At most max_percent, never all
        
        print(f"Constraints: min={min_expected}, max={max_expected} solvent accessible residues")
        
        # Add more residues if below minimum
        if len(self.solvent_accessible) < min_expected:
            print(f"Too few solvent-accessible residues detected ({len(self.solvent_accessible)}), "
                 f"adding more based on scores...")
            
            # Sort remaining residues by their exposure scores
            remaining = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.interacting_residues if r not in self.solvent_accessible],
                key=lambda x: x[1],  # Sort by score
                reverse=True  # Highest scores first
            )
            
            # Add only up to the minimum required
            for res_id, score in remaining:
                if len(self.solvent_accessible) >= min_expected:
                    break
                print(f"Adding {res_id} as solvent accessible based on ranking (score: {score:.2f})")
                self.solvent_accessible.add(res_id)
        
        # Remove residues if above maximum
        if len(self.solvent_accessible) > max_expected:
            print(f"Too many solvent-accessible residues detected ({len(self.solvent_accessible)}), "
                 f"removing lowest scoring ones...")
            
            # Sort current accessible residues by score, ascending
            to_evaluate = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.solvent_accessible],
                key=lambda x: x[1]  # Sort by score
            )
            
            # Remove lowest scoring residues until we're within limits
            residues_to_remove = len(self.solvent_accessible) - max_expected
            for i in range(residues_to_remove):
                if i < len(to_evaluate):
                    res_id, score = to_evaluate[i]
                    print(f"Removing {res_id} from solvent accessible (score: {score:.2f})")
                    self.solvent_accessible.remove(res_id)
        
        print(f"Final result: {len(self.solvent_accessible)} solvent-accessible residues out of {len(self.interacting_residues)} interacting residues")
        return self.solvent_accessible

    def calculate_dssp_solvent_accessibility(self, dssp_executable='dssp'):
        """
        Calculate solvent accessibility using DSSP.
        Requires DSSP executable to be installed and in PATH.
        
        Parameters:
        -----------
        dssp_executable : str
            Path to DSSP executable (default: 'dssp')
            
        Returns:
        --------
        set
            Set of (resname, resnum) tuples for solvent accessible residues
        """
        try:
            from Bio.PDB.DSSP import DSSP
            
            # Create a temporary PDB file for DSSP input
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_pdb:
                pdb_io = PDBIO()
                pdb_io.set_structure(self.structure)
                pdb_io.save(tmp_pdb.name)
                
                # Run DSSP
                dssp = DSSP(self.model, tmp_pdb.name, dssp=dssp_executable)
                
                # Clean up temporary file
                try:
                    os.unlink(tmp_pdb.name)
                except:
                    pass
                
                # Process DSSP results
                self.solvent_accessible = set()
                for (chain_id, res_id), dssp_data in dssp.property_dict.items():
                    resname = dssp_data[0]
                    resnum = res_id[1]
                    res_key = (resname, resnum)
                    
                    # Get relative solvent accessibility (0-1)
                    rsa = dssp_data[3]  # Relative accessibility
                    
                    # Consider residues with >25% accessibility as solvent accessible
                    if rsa > 0.25 and res_key in self.interacting_residues:
                        self.solvent_accessible.add(res_key)
                        
        except Exception as e:
            print(f"Warning: DSSP calculation failed. Falling back to geometric estimation. Error: {str(e)}")
            self.calculate_realistic_solvent_accessibility()
        
        return self.solvent_accessible
    
    def estimate_solvent_accessibility(self):
        """
        Estimate which residues might be solvent accessible using a simplified geometric approach.
        """
        # Try to use the realistic method first, fallback to simplified if it fails
        try:
            return self.calculate_realistic_solvent_accessibility()
        except Exception as e:
            print(f"Error in realistic solvent accessibility: {e}")
            print("Falling back to simplified estimation...")
            
            # For simplicity in the fallback, mark all residues as solvent accessible
            self.solvent_accessible = self.interacting_residues.copy()
            
            # But limit to at most 40% of residues
            if len(self.solvent_accessible) > 0:
                max_accessible = max(1, int(len(self.interacting_residues) * 0.4))
                if len(self.solvent_accessible) > max_accessible:
                    # Keep a random subset
                    self.solvent_accessible = set(list(self.solvent_accessible)[:max_accessible])
            
            return self.solvent_accessible
    
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
    
    def generate_interaction_report(self, output_file=None):
        """
        Generate a PLIP-like report for interactions.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the report will be saved
        
        Returns:
        --------
        str : Report text
        """
        if not output_file:
            base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
            output_file = f"{base_name}_interactions_report.txt"
        
        # Prepare ligand info
        if hasattr(self, 'ligand_residue') and self.ligand_residue:
            ligand_info = {
                'hetid': getattr(self.ligand_residue, 'resname', 'UNK'),
                'chain': 'X',
                'position': 0,
                'longname': getattr(self.ligand_residue, 'resname', 'Unknown'),
                'type': 'LIGAND',
                'interacting_chains': [],
                'interacting_res': []
            }
            
            # Try to get chain and position
            try:
                if hasattr(self.ligand_residue, 'parent') and self.ligand_residue.parent:
                    ligand_info['chain'] = self.ligand_residue.parent.id
                elif hasattr(self.ligand_residue, 'get_parent') and callable(self.ligand_residue.get_parent):
                    parent = self.ligand_residue.get_parent()
                    if parent and hasattr(parent, 'id'):
                        ligand_info['chain'] = parent.id
            except Exception as e:
                print(f"Error getting ligand chain: {str(e)}")
                
            try:
                if hasattr(self.ligand_residue, 'id') and isinstance(self.ligand_residue.id, tuple):
                    ligand_info['position'] = self.ligand_residue.id[1]
            except Exception as e:
                print(f"Error getting ligand position: {str(e)}")
        else:
            ligand_info = {
                'hetid': 'UNK', 
                'chain': 'X', 
                'position': 0, 
                'longname': 'Unknown Ligand', 
                'type': 'LIGAND',
                'interacting_chains': [], 
                'interacting_res': []
            }
        
        # Process interactions - tailored for our interactions structure
        processed_interactions = {}
        interacting_chains = set()
        interacting_res = set()
        
        for itype, interactions_list in self.interactions.items():
            processed_interactions[itype] = []
            
            for interaction in interactions_list:
                try:
                    # Extract basic interaction data
                    interaction_data = {
                        'distance': interaction['distance'],
                        'restype': interaction['protein_residue'].resname,
                        'resnr': interaction['protein_residue'].id[1],
                        'reschain': 'X'  # Default chain ID
                    }
                    
                    # Try to get chain
                    try:
                        if hasattr(interaction['protein_residue'], 'get_parent') and callable(interaction['protein_residue'].get_parent):
                            parent = interaction['protein_residue'].get_parent()
                            if parent and hasattr(parent, 'id'):
                                interaction_data['reschain'] = parent.id
                                interacting_chains.add(parent.id)
                    except Exception as e:
                        print(f"Error getting chain ID: {str(e)}")
                    
                    # Record interacting residue
                    res_id = f"{interaction_data['resnr']}{interaction_data['reschain']}"
                    interacting_res.add(res_id)
                    
                    processed_interactions[itype].append(interaction_data)
                    
                except Exception as e:
                    print(f"Warning: Error processing interaction: {str(e)}")
                    continue
        
        # Update ligand info with collected data
        ligand_info['interacting_chains'] = list(interacting_chains) if interacting_chains else []
        ligand_info['interacting_res'] = list(interacting_res) if interacting_res else []
        
        # Start building the report
        report = []
        report.append("=============================================================================")
        report.append(f"PandaMap-Color Interaction Report")
        report.append("=============================================================================")
        report.append("")
        report.append(f"Ligand: {ligand_info['hetid']}:{ligand_info['chain']}:{ligand_info['position']}")
        report.append(f"Name: {ligand_info['longname']}")
        report.append(f"Type: {ligand_info['type']}")
        report.append(f"File: {os.path.basename(self.pdb_file)}")
        report.append("\n------------------------------\n")
        
        # Summary statistics
        interacting_chains = ligand_info['interacting_chains']
        interacting_res_count = len(self.interacting_residues)
        report.append(f"Interacting Chains: {', '.join(interacting_chains) if interacting_chains else 'N/A'}")
        report.append(f"Interacting Residues: {interacting_res_count}")
        report.append(f"Solvent Accessible Residues: {len(self.solvent_accessible)}")
        
        # Summary of interactions
        interaction_types = {
            'hydrogen_bonds': "Hydrogen Bonds",
            'hydrophobic': "Hydrophobic Interactions",
            'pi_pi_stacking': "π-π Stacking",
            'carbon_pi': "Carbon-π Interactions",
            'donor_pi': "Donor-π Interactions",
            'amide_pi': "Amide-π Interactions",
            'halogen_bonds': "Halogen Bonds",
            'metal_coordination': "Metal Coordination",
            'ionic': "Ionic Interactions",
            'salt_bridge': "Salt Bridges",
            'covalent': "Covalent Bonds",
            'alkyl_pi': "Alkyl-π Interactions",
            'attractive_charge': "Attractive Charge",
            'pi_cation': "π-Cation Interactions",
            'repulsion': "Repulsion"
        }
        
        report.append("\n------------------------------\n")
        report.append("Interaction Summary:")
        for itype, label in interaction_types.items():
            count = len(self.interactions.get(itype, []))
            if count > 0:
                report.append(f"  {label}: {count}")
        
        # Detailed information about each interaction type
        for itype, label in interaction_types.items():
            if processed_interactions.get(itype, []):
                report.append("\n------------------------------\n")
                report.append(f"{label}:")
                for i, interaction in enumerate(processed_interactions[itype], 1):
                    # Format interaction details nicely
                    res_info = f"{interaction.get('restype', 'UNK')}{interaction.get('resnr', '?')}{interaction.get('reschain', '?')}"
                    dist_info = f"{interaction.get('distance', 0.0):.2f}Å"
                    lig_info = ligand_info.get('hetid', 'UNK')
                    
                    # Check if the residue is solvent accessible
                    res_key = (interaction.get('restype', 'UNK'), interaction.get('resnr', 0))
                    solvent_info = " (Solvent accessible)" if res_key in self.solvent_accessible else ""
                    
                    report.append(f"  {i}. {res_info}{solvent_info} -- {dist_info} -- {lig_info}")
        
        report.append("\n=============================================================================\n")
        
        # Write to file if specified
        if output_file:
            try:
                with open(output_file, 'w') as f:
                    f.write("\n".join(report))
                print(f"Report saved to {output_file}")
            except Exception as e:
                print(f"Error writing report to file: {str(e)}")
        
        return "\n".join(report)

    # visualization methods are now imported from visualization.py
   
    
    def run_analysis(self, output_file=None, use_dssp=False, generate_report=True, report_file=None):
        """
        Run the complete analysis pipeline with enhanced interaction detection.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved. If None, a default name will be generated.
        use_dssp : bool
            Whether to use DSSP for solvent accessibility (default: False)
        generate_report : bool
            Whether to generate a PLIP-like text report (default: True)
        report_file : str, optional
            Path where the report will be saved. If None but generate_report is True,
            a default path will be used.
            
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
            
            # Detect protein-ligand interactions with enhanced detection
            print("Detecting interactions...")
            self.detect_interactions()
            
            # Calculate solvent accessibility
            print("Calculating solvent accessibility...")
            if use_dssp:
                try:
                    self.calculate_dssp_solvent_accessibility()
                    
                    # Check if DSSP found reasonable number of solvent accessible residues
                    if len(self.solvent_accessible) > len(self.interacting_residues) * 0.5:
                        print("DSSP found too many solvent accessible residues, using realistic method")
                        self.calculate_realistic_solvent_accessibility()
                    elif len(self.solvent_accessible) < 2:
                        print("DSSP didn't find enough solvent accessible residues, using realistic method")
                        self.calculate_realistic_solvent_accessibility()
                except Exception as e:
                    print(f"Error using DSSP: {e}")
                    print("Falling back to realistic solvent accessibility calculation...")
                    self.calculate_realistic_solvent_accessibility()
            else:
                self.calculate_realistic_solvent_accessibility()
            
            # Generate visualization
            print("Generating visualization...")
            from .visualization import visualize
            result = visualize(self, output_file=output_file)
            
            # Generate text report if requested
            if generate_report:
                if report_file is None:
                    base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
                    report_file = f"{base_name}_interactions_report.txt"
                
                print("Generating interaction report...")
                try:
                    self.generate_interaction_report(output_file=report_file)
                except Exception as e:
                    print(f"Error generating report: {str(e)}")
            
            return result
            
        except Exception as e:
            error_msg = f"Error in analysis pipeline: {e}"
            print(f"✗ {error_msg}")
            return error_msg