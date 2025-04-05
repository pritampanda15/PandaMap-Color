#!/usr/bin/env python
"""
PandaMapColor: A Python package for visualizing protein-ligand 
interactions with customizable visual styling and design elements.

Visualization module for PandaMapColor.
Enhanced with improved interaction visualization from PandaMap.
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

# BioPython imports
from Bio.PDB import PDBParser, NeighborSearch, Selection

def visualize(self, output_file=None,
              figsize=(12, 12), dpi=300, title=None, 
              color_by_type=True, jitter=0.1,
              show_directionality=True):
    """
    Generate a 2D visualization of protein-ligand interactions
    with customizable styling.
    
    Parameters:
    -----------
    output_file : str, optional
        Path where the output image will be saved. If None, a default name will be generated.
    figsize : tuple
        Figure size in inches (width, height)
    dpi : int
        Resolution in dots per inch
    title : str, optional
        Title for the plot
    color_by_type : bool
        Whether to color residues by type (hydrophobic, polar, etc.)
    jitter : float
        Amount of random variation in position (0.0-1.0) for more natural look
    show_directionality : bool
        Whether to show interaction directionality with arrows (requires interaction_direction dict)
    """
    # If output_file is None, create a default filename
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
        output_file = f"{base_name}_interactions.png"
        print(f"No output file specified, using default: {output_file}")

    # Create the output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Set a clean style for plotting
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except:
        try:
            plt.style.use('seaborn-whitegrid')
        except:
            print("Warning: Could not set preferred plotting style. Using default style.")
            plt.style.use('default')
    
    # Debug check for solvent accessibility
    print("\n=== VISUALIZATION DEBUG ===")
    print(f"Total interacting residues: {len(self.interacting_residues)}")
    print(f"Solvent accessible residues: {len(self.solvent_accessible)}")
    
    if self.solvent_accessible:
        print("Solvent accessible residues:")
        for res_id in sorted(self.solvent_accessible):
            print(f"  - {res_id}")
    else:
        print("WARNING: No solvent accessible residues detected!")
    
    # Force reasonable solvent accessibility if all residues are marked or none are marked
    if len(self.solvent_accessible) == len(self.interacting_residues) and len(self.interacting_residues) > 0:
        print("WARNING: All residues are marked as solvent accessible!")
        print("This is likely incorrect. Removing some residues from solvent_accessible set...")
        
        # If all residues are marked, keep only about 40% 
        # Focus on residues that are typically exposed
        likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        keep_residues = []
        
        for res_id in self.interacting_residues:
            if res_id[0] in likely_exposed:
                keep_residues.append(res_id)
                
        # If we don't have enough exposed residues, add some hydrophobic ones near the surface
        max_to_keep = max(2, int(len(self.interacting_residues) * 0.4))
        if len(keep_residues) < max_to_keep:
            for res_id in self.interacting_residues:
                if res_id not in keep_residues and len(keep_residues) < max_to_keep:
                    keep_residues.append(res_id)
        
        # Replace the solvent_accessible set with our filtered version
        self.solvent_accessible = set(keep_residues[:max_to_keep])
        
        print(f"After filtering: {len(self.solvent_accessible)} solvent accessible residues")
        for res_id in sorted(self.solvent_accessible):
            print(f"  - {res_id}")
    
    # Prepare the title text
    if title:
        title_text = title
    else:
        title_text = f"Protein-Ligand Interactions: {os.path.basename(self.pdb_file)}"
    
    # Create a figure with dedicated space for title
    # Increase the figure height slightly to accommodate the title
    fig_height = figsize[1] + 1.0  # Add 1 inch for title
    fig = plt.figure(figsize=(figsize[0], fig_height), dpi=dpi)
    
    # Set figure background color
    fig.patch.set_facecolor('white')
    
    # Create a gridspec layout with 2 rows and 1 column
    # The top row (for title) will be smaller than the bottom row (for visualization)
    gs = plt.GridSpec(2, 1, height_ratios=[1, 10], figure=fig)
    
    # Create the title axes and main visualization axes
    title_ax = fig.add_subplot(gs[0])
    main_ax = fig.add_subplot(gs[1])
    
    # Configure the title axes
    title_ax.axis('off')  # Hide axes elements
    
    if self.colors.get('background_color', '#F8F8FF') == '#1A1A1A':  # If using dark mode
        # For dark mode, use normal weight (not bold) white text with subtle glow
        title_obj = title_ax.text(
            0.5, 0.5,  # Center position
            title_text,
            fontsize=18,
            fontweight='normal',  # Normal weight instead of bold for better readability
            ha='center',
            va='center',
            color='#FFFFFF'  # White text
        )
        
        # Add a subtle glow effect for dark mode
        title_obj.set_path_effects([
            path_effects.withStroke(linewidth=2.5, foreground='#333333'),
            path_effects.Normal()
        ])
    else:
        # For light mode, use bold black text
        title_obj = title_ax.text(
            0.5, 0.5,  # Center position
            title_text,
            fontsize=18,
            fontweight='bold',
            ha='center',
            va='center',
            color='black'
        )
        
        # Add standard shadow for light mode
        if self.use_enhanced_styling:
            title_obj.set_path_effects([
                path_effects.withStroke(linewidth=3, foreground='white'),
                path_effects.Normal()
            ])
    
    # Get background color from color scheme
    bg_color = self.colors.get('background_color', '#F8F8FF')
    
    # Create a subtle gradient background for the visualization area
    if self.colors.get('background_color', '#F8F8FF') == '#1A1A1A':  # If using dark mode
        # For dark mode, use full opacity
        background = Rectangle((-1000, -1000), 2000, 2000, 
                      facecolor='#1A1A1A', alpha=1.0, zorder=-1)
    else:
        # For light mode, use semi-transparent
        background = Rectangle((-1000, -1000), 2000, 2000, 
                      facecolor=bg_color, alpha=0.5, zorder=-1)
    main_ax.add_patch(background)
    
    # Add light blue background for ligand
    ligand_radius = 90
    ligand_pos = (0, 0)
    ligand_bg_color = self.colors.get('ligand_bg_color', '#ADD8E6')
    
    if self.use_enhanced_styling:
        # Inner glow effect for ligand area
        for r in np.linspace(ligand_radius, ligand_radius*0.5, 5):
            alpha = 0.1 * (1 - r/ligand_radius)
            glow = Circle(ligand_pos, r, facecolor=ligand_bg_color, 
                        edgecolor='none', alpha=alpha, zorder=0.5)
            main_ax.add_patch(glow)
    
    # Main ligand background
    ligand_circle = Circle(ligand_pos, ligand_radius, 
                         facecolor=ligand_bg_color, 
                         edgecolor='#87CEEB' if self.use_enhanced_styling else None, 
                         alpha=0.4, linewidth=1, zorder=1)
    main_ax.add_patch(ligand_circle)
    
    # Draw the ligand structure
    atom_positions = self.ligand_structure.draw_on_axes(main_ax, center=ligand_pos, radius=ligand_radius*0.8)
    
    # Place interacting residues in a circle around the ligand
    n_residues = len(self.interacting_residues)
    if n_residues == 0:
        print("Warning: No interacting residues detected.")
        n_residues = 1  # Avoid division by zero
        
    # Calculate positions for residues
    radius = 250  # Distance from center to residues
    residue_positions = {}
    
    # Get residue type colors
    residue_colors = self.colors.get('residue_colors', {})
    
    # Classify amino acids by type
    hydrophobic_aas = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'}
    polar_aas = {'SER', 'THR', 'ASN', 'GLN', 'CYS'}
    charged_pos_aas = {'LYS', 'ARG', 'HIS'}
    charged_neg_aas = {'ASP', 'GLU'}
    aromatic_aas = {'PHE', 'TYR', 'TRP'}
    
    # Arrange residues in a circle with slight randomization for natural look
    np.random.seed(42)  # For reproducibility
    
    # Dimensions for residue boxes
    rect_width, rect_height = 60, 40  # Residue box dimensions
    
    # For debugging
    print(f"Drawing {n_residues} residue nodes, {len(self.solvent_accessible)} with solvent accessibility")
    
    for i, res_id in enumerate(sorted(self.interacting_residues)):
        # Add slight randomness to angle and radius for more natural look
        angle = 2 * math.pi * i / n_residues
        angle_jitter = angle + np.random.uniform(-0.1, 0.1) * jitter
        radius_jitter = radius * (1 + np.random.uniform(-jitter, jitter))
        
        x = radius_jitter * math.cos(angle_jitter)
        y = radius_jitter * math.sin(angle_jitter)
        residue_positions[res_id] = (x, y)
        
        # Draw solvent accessibility highlight
        if res_id in self.solvent_accessible:
            print(f"Drawing solvent accessibility circle for {res_id}")
            solvent_color = self.colors.get('solvent_color', '#ADD8E6')
            
            if self.use_enhanced_styling:
                # Create a glow effect for solvent accessibility
                for r in np.linspace(40, 20, 3):
                    alpha = 0.1 * (1 - r/40)
                    glow = Circle((x, y), r, facecolor=solvent_color, 
                               edgecolor=None, alpha=alpha, zorder=0.8)
                    main_ax.add_patch(glow)
            
            solvent_circle = Circle((x, y), 40, facecolor=solvent_color, 
                                  edgecolor='#87CEEB', alpha=0.5, zorder=1)
            main_ax.add_patch(solvent_circle)
        else:
            print(f"Residue {res_id} is NOT solvent accessible - no circle")
        
        # Determine residue type color
        resname, resnum = res_id
        if color_by_type:
            if resname in hydrophobic_aas:
                res_color = residue_colors.get('hydrophobic', '#FFD700')
                edge_color = '#B8860B'  # Dark goldenrod
            elif resname in polar_aas:
                res_color = residue_colors.get('polar', '#00BFFF')
                edge_color = '#0000CD'  # Medium blue
            elif resname in charged_pos_aas:
                res_color = residue_colors.get('charged_pos', '#FF6347')
                edge_color = '#8B0000'  # Dark red
            elif resname in charged_neg_aas:
                res_color = residue_colors.get('charged_neg', '#32CD32')
                edge_color = '#006400'  # Dark green
            elif resname in aromatic_aas:
                res_color = residue_colors.get('aromatic', '#DA70D6')
                edge_color = '#8B008B'  # Dark magenta
            else:
                res_color = residue_colors.get('default', 'white')
                edge_color = 'black'
        else:
            res_color = residue_colors.get('default', 'white')
            edge_color = 'black'
        
        if self.use_enhanced_styling:
            # Create a subtle shadow for 3D effect
            shadow_offset = 2
            shadow = self.create_rounded_rectangle(
                (x-30+shadow_offset, y-20+shadow_offset), 60, 40, radius=10,
                facecolor='gray', edgecolor=None, alpha=0.2, zorder=1.9
            )
            main_ax.add_patch(shadow)
        
        # Draw residue node as rounded rectangle
        if self.use_enhanced_styling:
            residue_box = self.create_rounded_rectangle(
                (x-30, y-20), 60, 40, radius=10,
                facecolor=res_color, edgecolor=edge_color, linewidth=1.5,
                alpha=0.9, zorder=2
            )
            
            # Add a subtle inner highlight for a 3D effect
            highlight = self.create_rounded_rectangle(
                (x-27, y-17), 54, 34, radius=8,
                facecolor='white', edgecolor=None, alpha=0.2, zorder=2.1
            )
            main_ax.add_patch(highlight)
        else:
            # Simpler rectangular node
            residue_box = Rectangle((x-30, y-20), 60, 40,
                                  facecolor=res_color, edgecolor=edge_color, linewidth=1.5,
                                  alpha=0.9, zorder=2)
        
        main_ax.add_patch(residue_box)
        
        # Add residue label (NAME NUMBER)
        resname, resnum = res_id
        label = f"{resname} {resnum}"
        text = main_ax.text(x, y, label, ha='center', va='center',
                          fontsize=11, fontweight='bold', zorder=3)
        
        # Add text shadow for better readability
        if self.use_enhanced_styling:
            text.set_path_effects([
                path_effects.withStroke(linewidth=3, foreground='white'),
                path_effects.Normal()
            ])
        else:
            text.set_path_effects([
                path_effects.withStroke(linewidth=2, foreground='white')
            ])
    
    # Get interaction styles from color scheme
    interaction_styles = self.colors.get('interaction_styles', {})
    
    # Add additional interaction styles from PandaMap if they don't exist
    additional_styles = {
        'ionic': {'color': '#FF4500', 'linestyle': '-', 'linewidth': 1.5,
                'marker_text': 'I', 'marker_color': '#FF4500', 'marker_bg': '#FFE4E1',
                'name': 'Ionic', 'glow': False},
        'halogen_bonds': {'color': '#00CED1', 'linestyle': '-', 'linewidth': 1.5,
                       'marker_text': 'X', 'marker_color': '#00CED1', 'marker_bg': '#E0FFFF',
                       'name': 'Halogen Bond', 'glow': False},
        'cation_pi': {'color': '#FF00FF', 'linestyle': '--', 'linewidth': 1.5,
                    'marker_text': 'C+π', 'marker_color': '#FF00FF', 'marker_bg': 'white',
                    'name': 'Cation-Pi', 'glow': False},
        'metal_coordination': {'color': '#FFD700', 'linestyle': '-', 'linewidth': 1.5,
                            'marker_text': 'M', 'marker_color': '#FFD700', 'marker_bg': '#FFFACD',
                            'name': 'Metal Coordination', 'glow': False},
        'salt_bridge': {'color': '#FF6347', 'linestyle': '-', 'linewidth': 1.5,
                      'marker_text': 'S', 'marker_color': '#FF6347', 'marker_bg': '#FFEFD5',
                      'name': 'Salt Bridge', 'glow': False},
        'covalent': {'color': '#000000', 'linestyle': '-', 'linewidth': 2.0,
                   'marker_text': 'COV', 'marker_color': '#000000', 'marker_bg': '#FFFFFF',
                   'name': 'Covalent Bond', 'glow': False},
        'alkyl_pi': {'color': '#4682B4', 'linestyle': '--', 'linewidth': 1.5,
                   'marker_text': 'A-π', 'marker_color': '#4682B4', 'marker_bg': 'white',
                   'name': 'Alkyl-Pi', 'glow': False},
        'attractive_charge': {'color': '#1E90FF', 'linestyle': '-', 'linewidth': 1.5,
                           'marker_text': 'A+', 'marker_color': '#1E90FF', 'marker_bg': '#E6E6FA',
                           'name': 'Attractive Charge', 'glow': False},
        'pi_cation': {'color': '#FF00FF', 'linestyle': '--', 'linewidth': 1.5,
                    'marker_text': 'π-C+', 'marker_color': '#FF00FF', 'marker_bg': 'white',
                    'name': 'Pi-Cation', 'glow': False},
        'repulsion': {'color': '#DC143C', 'linestyle': '-', 'linewidth': 1.5,
                    'marker_text': 'R', 'marker_color': '#DC143C', 'marker_bg': '#FFC0CB',
                    'name': 'Repulsion', 'glow': False}
    }
    
    # Add new interaction styles if not already present
    for int_type, style in additional_styles.items():
        if int_type not in interaction_styles:
            interaction_styles[int_type] = style
    
    # Update marker text to match reference image
    if 'carbon_pi' in interaction_styles:
        interaction_styles['carbon_pi']['marker_text'] = 'C-π'
    if 'pi_pi_stacking' in interaction_styles:
        interaction_styles['pi_pi_stacking']['marker_text'] = 'π-π'
    
    # Function to find box edge intersection
    def find_box_edge(box_center, target_point, width, height):
        """Find where a line from box center to target point intersects the box edge"""
        dx = target_point[0] - box_center[0]
        dy = target_point[1] - box_center[1]
        angle = math.atan2(dy, dx)
        
        half_width = width/2
        half_height = height/2
        
        if abs(dx) > abs(dy):
            x_intersect = box_center[0] + (half_width if dx > 0 else -half_width)
            y_intersect = box_center[1] + (x_intersect - box_center[0]) * dy/dx
            if abs(y_intersect - box_center[1]) > half_height:
                y_intersect = box_center[1] + (half_height if dy > 0 else -half_height)
                x_intersect = box_center[0] + (y_intersect - box_center[1]) * dx/dy
        else:
            y_intersect = box_center[1] + (half_height if dy > 0 else -half_height)
            x_intersect = box_center[0] + (y_intersect - box_center[1]) * dx/dy
            if abs(x_intersect - box_center[0]) > half_width:
                x_intersect = box_center[0] + (half_width if dx > 0 else -half_width)
                y_intersect = box_center[1] + (x_intersect - box_center[0]) * dy/dx
                
        return (x_intersect, y_intersect)
    
    # Store interaction lines for marker placement
    interaction_lines = []
    
    # First pass: group interactions by residue to better handle overlaps
    interactions_by_residue = defaultdict(list)
    for interaction_type, interactions in self.interactions.items():
        if interaction_type not in interaction_styles:
            continue
        
        for interaction in interactions:
            res = interaction['protein_residue']
            res_id = (res.resname, res.id[1])
            if res_id in residue_positions:
                interactions_by_residue[res_id].append({
                    'interaction_type': interaction_type,
                    'interaction': interaction
                })
    
    # Draw interaction lines for each residue
    for res_id, grouped_interactions in interactions_by_residue.items():
        res_pos = residue_positions[res_id]
        
        # If multiple interactions for same residue, we need to space them
        if len(grouped_interactions) > 1:
            # Calculate base angle from residue to ligand center
            base_angle = math.atan2(res_pos[1] - ligand_pos[1], res_pos[0] - ligand_pos[0])
            
            # Calculate marker spacing in radians based on number of interactions
            angle_spacing = min(0.2, 0.6 / len(grouped_interactions))
            
            # Start with an offset to center the markers
            start_offset = -(len(grouped_interactions) - 1) * angle_spacing / 2
        else:
            base_angle = 0
            angle_spacing = 0
            start_offset = 0
        
        # Process each interaction for this residue
        for i, interaction_data in enumerate(grouped_interactions):
            interaction_type = interaction_data['interaction_type']
            interaction = interaction_data['interaction']
            lig_atom = interaction['ligand_atom']
            style = interaction_styles[interaction_type]
            use_glow = style.get('glow', False) and self.use_enhanced_styling
            
            # Try to use actual ligand atom position if available
            if lig_atom.get_id() in atom_positions:
                lig_pos = atom_positions[lig_atom.get_id()]
            else:
                # Determine the point on ligand circle edge as fallback
                angle = base_angle
                lig_edge_x = ligand_pos[0] + ligand_radius * math.cos(angle)
                lig_edge_y = ligand_pos[1] + ligand_radius * math.sin(angle)
                lig_pos = (lig_edge_x, lig_edge_y)
            
            # Find box edge intersection for better line placement
            box_edge_pos = find_box_edge(res_pos, lig_pos, rect_width, rect_height)
            
            # Calculate angle offset for this interaction
            current_angle_offset = start_offset + i * angle_spacing
            
            # Calculate curvature
            dx = res_pos[0] - lig_pos[0]
            dy = res_pos[1] - lig_pos[1]
            distance = math.sqrt(dx*dx + dy*dy)
            curvature = 0.08 * (200 / max(distance, 100))
            
            # Add a subtle glow effect to important interactions in standard mode
            if use_glow:
                # Draw a thicker, semi-transparent line underneath
                glow = FancyArrowPatch(
                    lig_pos, res_pos,
                    connectionstyle=f"arc3,rad={curvature}",
                    color=style['color'],
                    linestyle=style['linestyle'],
                    linewidth=style['linewidth'] * 2.5,  # Thicker for glow
                    arrowstyle='-',
                    alpha=0.2,  # Semi-transparent
                    zorder=3.8
                )
                main_ax.add_patch(glow)
            
            # Store line parameters for marker placement
            line_params = {
                'start_pos': box_edge_pos,
                'end_pos': lig_pos,
                'curvature': curvature,
                'style': style,
                'interaction_type': interaction_type,
                'res_id': res_id,
                'key': f"{interaction_type}_{res_id[0]}_{res_id[1]}",
                'distance': distance
            }
            interaction_lines.append(line_params)
            
            # Check if we should show directionality
            if show_directionality and hasattr(self, 'interaction_direction'):
                # Get directionality info
                interaction_key = (res_id, interaction_type)
                direction = self.interaction_direction.get(interaction_key, 'bidirectional')
                
                # Draw appropriate arrow based on direction
                if direction == 'protein_to_ligand':
                    # Arrow from protein to ligand
                    arrow = FancyArrowPatch(
                        box_edge_pos, lig_pos,
                        connectionstyle=f"arc3,rad={curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='-|>',
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    main_ax.add_patch(arrow)
                    
                elif direction == 'ligand_to_protein':
                    # Arrow from ligand to protein
                    arrow = FancyArrowPatch(
                        lig_pos, box_edge_pos,
                        connectionstyle=f"arc3,rad={-curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='-|>',
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    main_ax.add_patch(arrow)
                    
                else:  # bidirectional - single line with arrows on both ends
                    # Draw one bidirectional arrow
                    arrow = FancyArrowPatch(
                        box_edge_pos, lig_pos,
                        connectionstyle=f"arc3,rad={curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='<|-|>',  # Arrows on both ends
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    main_ax.add_patch(arrow)
            else:
                # Standard non-directional line
                line = FancyArrowPatch(
                    box_edge_pos, lig_pos,
                    connectionstyle=f"arc3,rad={curvature}",
                    color=style['color'],
                    linestyle=style['linestyle'],
                    linewidth=style['linewidth'],
                    arrowstyle='-',
                    alpha=0.85 if self.use_enhanced_styling else 0.7,
                    zorder=4,
                    capstyle='round' if self.use_enhanced_styling else 'butt',
                    joinstyle='round' if self.use_enhanced_styling else 'miter'
                )
                main_ax.add_patch(line)
    
    # Calculate marker positions along the interaction lines
    marker_positions = {}
    
    # Sort interactions by type for consistent placement
    # Place hydrogen bonds first, then pi interactions, then hydrophobic
    type_order = {
        'hydrogen_bonds': 0, 'ionic': 1, 'salt_bridge': 2, 'halogen_bonds': 3,
        'metal_coordination': 4, 'pi_pi_stacking': 5, 'cation_pi': 6, 
        'carbon_pi': 7, 'donor_pi': 8, 'amide_pi': 9, 'hydrophobic': 10,
        'alkyl_pi': 11, 'attractive_charge': 12, 'pi_cation': 13, 'repulsion': 14
    }
    
    sorted_lines = sorted(interaction_lines,
                        key=lambda x: (type_order.get(x['interaction_type'], 999), x['distance']))
    
    # Second pass: Place markers along paths
    for line_params in sorted_lines:
        start_pos = line_params['start_pos']
        end_pos = line_params['end_pos']
        curvature = line_params['curvature']
        style = line_params['style']
        key = line_params['key']
        res_id = line_params['res_id']
        interaction_type = line_params['interaction_type']
        
        # Get directionality for adjusting marker position
        direction = 'bidirectional'
        if hasattr(self, 'interaction_direction'):
            direction = self.interaction_direction.get((res_id, interaction_type), 'bidirectional')
        
        # Calculate points along the curved path
        path_points = []
        steps = 20
        for i in range(steps + 1):
            t = i / steps  # Parameter along the curve (0 to 1)
            
            # For bidirectional, use the protein->ligand curve for marker
            # For directional, offset slightly based on direction
            curve_adjust = curvature
            if direction == 'ligand_to_protein':
                curve_adjust = -curvature
            
            # Quadratic Bezier curve formula for approximating arc
            control_x = (start_pos[0] + end_pos[0])/2 + curve_adjust * (end_pos[1] - start_pos[1]) * 2
            control_y = (start_pos[1] + end_pos[1])/2 - curve_adjust * (end_pos[0] - start_pos[0]) * 2
            
            # Calculate point at parameter t
            x = (1-t)*(1-t)*start_pos[0] + 2*(1-t)*t*control_x + t*t*end_pos[0]
            y = (1-t)*(1-t)*start_pos[1] + 2*(1-t)*t*control_y + t*t*end_pos[1]
            
            path_points.append((x, y))
        
        # Find best marker position
        best_position = None
        best_score = float('-inf')
        
        # Try different positions along the path to find the best placement
        t_values = [0.5, 0.45, 0.55, 0.4, 0.6, 0.35, 0.65, 0.3, 0.7, 0.25, 0.75]
        
        # Special case for C-π, use a specific t value that's known to work well
        if interaction_type == 'carbon_pi':
            t_values = [0.42, 0.4, 0.45, 0.38]  # Favors positions closer to the ligand
        
        for t in t_values:
            idx = int(t * steps)
            if idx >= len(path_points):
                idx = len(path_points) - 1
            pos = path_points[idx]
            
            # Calculate distance to existing markers
            if marker_positions:  # Only if there are existing markers
                min_dist = min(math.sqrt((pos[0]-p[0])**2 + (pos[1]-p[1])**2) 
                            for p in marker_positions.values())
            else:
                min_dist = float('inf')
            
            text_len = len(style['marker_text'])
            min_req_dist = 25 + text_len * 2
            score = min(min_dist / min_req_dist, 2.0) + (1.0 - abs(t - 0.5))
            
            if score > best_score:
                best_score = score
                best_position = pos
        
        if best_position is None:
            best_position = path_points[len(path_points)//2]
        
        marker_positions[key] = best_position
        x, y = best_position
        
        # Draw marker shape
        marker_radius = 9 + (len(style['marker_text']) - 1) * 1.5
        if 'pi' in interaction_type:
            angles = np.linspace(0, 2*np.pi, 7)[:-1]
            vertices = [(x + marker_radius * math.cos(a), y + marker_radius * math.sin(a)) 
                      for a in angles]
            marker = Polygon(vertices, closed=True, facecolor=style.get('marker_bg', 'white'),
                           edgecolor=style['color'], linewidth=1.5, zorder=5)
            main_ax.add_patch(marker)
        else:
            marker = Circle((x, y), marker_radius, facecolor=style.get('marker_bg', 'white'),
                          edgecolor=style['color'], linewidth=1.5, zorder=5)
            main_ax.add_patch(marker)
        
        # Add marker text
        font_size = 9 if len(style['marker_text']) <= 1 else max(7, 9 - (len(style['marker_text']) - 1) * 0.8)
        text = main_ax.text(x, y, style['marker_text'], ha='center', va='center',
                          fontsize=font_size, color=style['color'], fontweight='bold', zorder=6)
        text.set_path_effects([path_effects.withStroke(linewidth=1, foreground='white')])
    
    # Create legend for interaction types
    legend_title = "Interacting structural groups"
    legend_elements = []

    # Add residue example for legend
    residue_patch = Rectangle((0, 0), 1, 1, facecolor='white', 
                            edgecolor='black', label='Interacting structural groups')
    legend_elements.append(residue_patch)
    
    # Interaction type markers for legend
    for int_type, style in sorted(interaction_styles.items(), 
                                key=lambda x: type_order.get(x[0], 999)):
        # Only include interaction types that are present
        if int_type in self.interactions and self.interactions[int_type]:
            if int_type == 'hydrogen_bonds':
                # Special handling for H-bonds to match reference
                line = Line2D([0], [0], color=style['color'],
                            linestyle=style['linestyle'], linewidth=style['linewidth'],
                            marker='o', markerfacecolor=style.get('marker_bg', 'white'), 
                            markeredgecolor=style['color'],
                            markersize=8, label=style['name'])
            elif 'pi' in int_type:
                # Hexagonal markers for pi interactions
                line = Line2D([0], [0], color=style['color'],
                            linestyle=style['linestyle'], linewidth=style['linewidth'],
                            marker='h', markerfacecolor=style.get('marker_bg', 'white'), 
                            markeredgecolor=style['color'],
                            markersize=8, label=style['name'])
            else:
                # Circular markers for other interactions
                line = Line2D([0], [0], color=style['color'],
                            linestyle=style['linestyle'], linewidth=style['linewidth'],
                            marker='o', markerfacecolor=style.get('marker_bg', 'white'), 
                            markeredgecolor=style['color'],
                            markersize=8, label=style['name'])
            legend_elements.append(line)
    
    # Add directionality to legend if showing directionality
    if show_directionality and hasattr(self, 'interaction_direction'):
        legend_elements.append(
            Line2D([0], [0], color='black', linestyle='-', marker='>',
                markerfacecolor='black', markersize=8, 
                label='Unidirectional interaction')
        )
        # For bidirectional arrow, use a simpler approach with a diamond symbol
        bidirectional = Line2D([0], [0], color='black', linestyle='-',
                            marker='d', markerfacecolor='black',
                            markersize=8, label='Bidirectional interaction')
        legend_elements.append(bidirectional)
    
    # Add solvent accessibility indicator
    if self.solvent_accessible:
        solvent_color = self.colors.get('solvent_color', '#ADD8E6')
        solvent_patch = Rectangle((0, 0), 1, 1, facecolor=solvent_color, 
                                alpha=0.5, edgecolor='#87CEEB', label='Solvent accessible')
        legend_elements.append(solvent_patch)
    
    # Create an enhanced legend box in top right corner
    legend_font_size = 10 if self.use_enhanced_styling else 9
    title_font_size = 11 if self.use_enhanced_styling else 10
    
    if self.colors.get('background_color', '#F8F8FF') == '#1A1A1A':  # If using dark mode
        # Dark mode legend
        legend = main_ax.legend(
            handles=legend_elements,
            title=legend_title,
            loc='upper right',
            frameon=True,
            framealpha=0.95,  # More opaque for better contrast
            fontsize=legend_font_size,
            title_fontsize=title_font_size,
            facecolor='#2A2A2A',  # Slightly lighter than background
            edgecolor='#505050'  # Visible border
        )
    
        # Make legend text white for dark mode
        for text in legend.get_texts():
            text.set_color('black')  # Use white text for dark mode
    
        # Make legend title white and bold for dark mode
        title = legend.get_title()
        title.set_color('black')  # Use white title for dark mode
        
    else:
        # Light mode legend
        legend = main_ax.legend(
            handles=legend_elements,
            title=legend_title,
            loc='upper right',
            frameon=True,
            framealpha=0.85 if self.use_enhanced_styling else 0.7,
            fontsize=legend_font_size,
            title_fontsize=title_font_size,
            facecolor='white',
            edgecolor='gray'
        )
    
    # Standard legend styling for light mode
    if self.use_enhanced_styling:
        # Add a shadow effect to the legend
        frame = legend.get_frame()
        frame.set_linewidth(1.0)
        frame.set_facecolor('white')
        
        # Make the legend title bold
        title = legend.get_title()
        title.set_fontweight('bold')
    
    # Set plot limits and appearance
    max_coord = radius + 100
    main_ax.set_xlim(-max_coord, max_coord)
    main_ax.set_ylim(-max_coord, max_coord)
    main_ax.set_aspect('equal')
    main_ax.axis('off')
    
    # Adjust layout to ensure proper spacing
    plt.tight_layout()

    # Save the figure with the properly rendered title
    if self.colors.get('background_color', '#F8F8FF') == '#1A1A1A':  # If using dark mode
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', 
                  facecolor='#1A1A1A', edgecolor='none',
                  transparent=False, pad_inches=0.3)
    else:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', 
                  facecolor='white', edgecolor='none',
                  transparent=False, pad_inches=0.3)
    plt.close()

    print(f"Successfully saved visualization to: {output_file}")
    return output_file
    
def run_analysis(self, output_file=None, use_dssp=False, generate_report=True, report_file=None):
    """
    Run the complete analysis pipeline.
    
    Parameters:
    -----------
    output_file : str, optional
        Path where the output image will be saved. If None, a default name will be generated.
    use_dssp : bool
        Whether to use DSSP for solvent accessibility (default: False)
    generate_report : bool
        Whether to generate a text report (default: True)
    report_file : str, optional
        Path where the report will be saved
        
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
        
        # Detect protein-ligand interactions with enhanced detection if available
        print("Detecting interactions...")
        # Check if the improved interaction detection is available
        if hasattr(self, 'detect_interactions'):
            self.detect_interactions()
        
        # Calculate solvent accessibility
        print("Calculating solvent accessibility...")
        if use_dssp and hasattr(self, 'calculate_dssp_solvent_accessibility'):
            try:
                self.calculate_dssp_solvent_accessibility()
                
                # Check if DSSP found reasonable number of solvent accessible residues
                if len(self.solvent_accessible) > len(self.interacting_residues) * 0.5:
                    print("DSSP found too many solvent accessible residues, using realistic method")
                    
                    if hasattr(self, 'calculate_realistic_solvent_accessibility'):
                        self.calculate_realistic_solvent_accessibility()
                    else:
                        self.estimate_solvent_accessibility()
                        
                elif len(self.solvent_accessible) < 2 and len(self.interacting_residues) > 2:
                    print("DSSP didn't find enough solvent accessible residues, using realistic method")
                    
                    if hasattr(self, 'calculate_realistic_solvent_accessibility'):
                        self.calculate_realistic_solvent_accessibility()
                    else:
                        self.estimate_solvent_accessibility()
            except Exception as e:
                print(f"Error using DSSP: {e}")
                print("Falling back to alternative solvent accessibility calculation...")
                
                if hasattr(self, 'calculate_realistic_solvent_accessibility'):
                    self.calculate_realistic_solvent_accessibility()
                else:
                    self.estimate_solvent_accessibility()
        elif hasattr(self, 'calculate_realistic_solvent_accessibility'):
            self.calculate_realistic_solvent_accessibility()
        else:
            self.estimate_solvent_accessibility()
        
        # Generate visualization with directionality if available
        show_directionality = hasattr(self, 'interaction_direction')
        print("Generating visualization...")
        result = visualize(self, output_file=output_file, show_directionality=show_directionality)
        
        # Generate text report if requested
        if generate_report and hasattr(self, 'generate_interaction_report'):
            if report_file is None:
                base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
                report_file = f"{base_name}_interactions_report.txt"
            
            print("Generating interaction report...")
            try:
                self.generate_interaction_report(output_file=report_file)
                print(f"Report saved to {report_file}")
            except Exception as e:
                print(f"Error generating report: {str(e)}")
        
        return result
        
    except Exception as e:
        error_msg = f"Error in analysis pipeline: {e}"
        print(f"✗ {error_msg}")
        return error_msg