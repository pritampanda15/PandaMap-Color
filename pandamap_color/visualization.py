#!/usr/bin/env python
"""
PandaMapColor: A Python package for visualizing protein-ligand 
interactions with customizable visual styling and design elements.

Visualization module for PandaMapColor.
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
                  color_by_type=True, jitter=0.1):
    
            
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
            solvent_color = self.colors.get('solvent_color', '#ADD8E6')
            
            if self.use_enhanced_styling:
                # Create a glow effect for solvent accessibility
                for r in np.linspace(40, 20, 3):
                    alpha = 0.1 * (1 - r/40)
                    glow = Circle((x, y), r, facecolor=solvent_color, 
                                edgecolor=None, alpha=alpha, zorder=0.8)
                    main_ax.add_patch(glow)
            
            solvent_circle = Circle((x, y), 40, facecolor=solvent_color, 
                                    edgecolor='none', alpha=0.3, zorder=1)
            main_ax.add_patch(solvent_circle)
        
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
    
    # Update marker text to match reference image
    if 'carbon_pi' in interaction_styles:
        interaction_styles['carbon_pi']['marker_text'] = 'C-π'
    if 'pi_pi_stacking' in interaction_styles:
        interaction_styles['pi_pi_stacking']['marker_text'] = 'π-π'
    
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
    
    # Dictionary to track all interaction lines and their parameters
    interaction_lines = {}
    
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
            
            # Calculate angle offset for this interaction
            current_angle_offset = start_offset + i * angle_spacing
            
            # Draw curved line from ligand to residue
            # Adjust curvature based on position in the group
            curvature = 0.1
            if len(grouped_interactions) > 1:
                curvature = 0.1 + current_angle_offset  # Vary curvature to separate lines
            
            # Add a subtle glow effect to important interactions
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
            
            # Draw curved line from ligand to residue
            line = FancyArrowPatch(
                lig_pos, res_pos,
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
            
            # Store the line parameters for marker placement - IMPORTANT: store full interaction_type
            key = f"{interaction_type}_{res_id[0]}_{res_id[1]}"
            interaction_lines[key] = {
                'lig_pos': lig_pos,
                'res_pos': res_pos,
                'curvature': curvature,
                'style': style,
                'interaction_type': interaction_type  # Store the FULL interaction_type
            }
    
    # Calculate marker positions along the interaction lines
    marker_positions = {}
    
    # Sort interactions by type for consistent placement
    # Place hydrogen bonds first, then pi interactions, then hydrophobic
    type_order = {'hydrogen_bonds': 0, 'carbon_pi': 1, 'pi_pi_stacking': 2, 
                    'donor_pi': 3, 'amide_pi': 4, 'hydrophobic': 5}
    
    sorted_lines = sorted(interaction_lines.items(), 
                            key=lambda x: type_order.get(x[1]['interaction_type'], 999))
    
    # Second pass: Place markers along paths
    for key, line_params in sorted_lines:
        # Get parameters from the stored line data
        lig_pos = line_params['lig_pos']
        res_pos = line_params['res_pos']
        curvature = line_params['curvature']
        style = line_params['style']
        interaction_type = line_params['interaction_type']  # Use the stored full interaction type
        
        # Calculate the path of the curved line
        path_points = []
        steps = 20
        for i in range(steps + 1):
            t = i / steps  # Parameter along the curve (0 to 1)
            
            # Quadratic Bezier curve formula for approximating arc
            control_x = (lig_pos[0] + res_pos[0])/2 + curvature * (res_pos[1] - lig_pos[1]) * 2
            control_y = (lig_pos[1] + res_pos[1])/2 - curvature * (res_pos[0] - lig_pos[0]) * 2
            
            # Calculate point at parameter t
            x = (1-t)*(1-t)*lig_pos[0] + 2*(1-t)*t*control_x + t*t*res_pos[0]
            y = (1-t)*(1-t)*lig_pos[1] + 2*(1-t)*t*control_y + t*t*res_pos[1]
            
            path_points.append((x, y))
        
        # Try different positions along the path until finding one that doesn't overlap
        # Start with middle and work outward
        t_values = [0.5, 0.45, 0.55, 0.4, 0.6, 0.35, 0.65, 0.3, 0.7, 0.25, 0.75]
        
        # Special case for C-π, use a specific t value that's known to work well
        if interaction_type == 'carbon_pi':
            t_values = [0.42, 0.4, 0.45, 0.38]  # Favors positions closer to the ligand
        
        marker_placed = False
        best_position = None
        best_score = float('-inf')
        
        for t in t_values:
            idx = int(t * steps)
            if idx >= len(path_points):
                idx = len(path_points) - 1
            pos = path_points[idx]
            
            # Calculate distance to existing markers
            min_dist_to_markers = float('inf')
            for other_pos in marker_positions.values():
                dist = math.sqrt((pos[0] - other_pos[0])**2 + (pos[1] - other_pos[1])**2)
                min_dist_to_markers = min(min_dist_to_markers, dist)
                
            # Calculate distance to the line
            min_dist_to_line = float('inf')
            for i in range(len(path_points)-1):
                x1, y1 = path_points[i]
                x2, y2 = path_points[i+1]
                
                # Distance from point to line segment
                px, py = pos
                
                # Calculate projection
                line_length_sq = (x2-x1)**2 + (y2-y1)**2
                if line_length_sq == 0:
                    # Point 1 and 2 are the same
                    dist_to_segment = math.sqrt((px-x1)**2 + (py-y1)**2)
                else:
                    # Calculate projection parameter
                    t_proj = max(0, min(1, ((px-x1)*(x2-x1) + (py-y1)*(y2-y1)) / line_length_sq))
                    
                    # Find closest point on segment
                    closest_x = x1 + t_proj * (x2-x1)
                    closest_y = y1 + t_proj * (y2-y1)
                    
                    # Distance to that point
                    dist_to_segment = math.sqrt((px-closest_x)**2 + (py-closest_y)**2)
                
                min_dist_to_line = min(min_dist_to_line, dist_to_segment)
            
            # Adjust min distance based on marker text length
            text_length = len(style['marker_text'])
            min_marker_distance = 28 + text_length * 2.5  # Slightly reduced from earlier
            
            # Calculate a score that balances:
            # 1. Distance from other markers (higher is better)
            # 2. Proximity to the line (lower is better)
            # 3. Closeness to middle of the line (t near 0.5 is better)
            
            overlap_score = min(min_dist_to_markers / min_marker_distance, 2.0)  # Cap at 2.0
            line_proximity_score = 10.0 / (min_dist_to_line + 1.0)  # Higher for closer to line
            middle_preference = 1.0 - abs(t - 0.5) * 0.8  # Higher for positions near middle
            
            # If there's a major overlap, heavily penalize
            if min_dist_to_markers < min_marker_distance * 0.7:
                overlap_score = overlap_score * 0.2
            
            # Calculate total score - weighted sum
            total_score = (
                overlap_score * 1.0 +         # Avoid overlaps
                line_proximity_score * 3.0 +  # Strongly prefer positions near the line
                middle_preference * 0.5       # Slight preference for middle
            )
            
            # Update best position if score is higher
            if total_score > best_score:
                best_score = total_score
                best_position = pos
        
        # If no non-overlapping position, use the one with least overlap
        if not marker_placed and best_position:
            pass  # Already stored in best_position
        
        # If still no position, use midpoint as fallback
        if not best_position:
            best_position = path_points[int(len(path_points)/2)]
        
        # Store the final position
        marker_positions[key] = best_position
    
    # Now draw all the markers at the calculated positions
    for key, position in marker_positions.items():
        # Extract the interaction type from the key - but don't rely on splitting
        # Since the key format is "{interaction_type}_{res_id[0]}_{res_id[1]}",
        # we need to get the original interaction_type
        interaction_type = interaction_lines[key]['interaction_type']
        
        # Get the style for this interaction type
        style = interaction_styles[interaction_type]
        
        # Set the position
        midpoint_x, midpoint_y = position
        
        # Draw marker based on interaction type
        n_sides = 6 if 'pi' in interaction_type else 0  # Hexagon for pi interactions
        
        # Adjust marker radius based on text length
        text_length = len(style['marker_text'])
        if 'pi' in interaction_type and text_length > 1:
            marker_radius = 14  # Larger for multi-character pi interactions
        else:
            marker_radius = 12  # Standard for single-character markers
        
        if n_sides > 0:
            # Draw a hexagon for pi interactions
            angles = np.linspace(0, 2*np.pi, n_sides+1)[:-1]
            hex_vertices = np.array([
                [midpoint_x + marker_radius * np.cos(angle),
                    midpoint_y + marker_radius * np.sin(angle)]
                for angle in angles
            ])
            hex_patch = Polygon(
                hex_vertices, 
                closed=True,
                facecolor=style.get('marker_bg', 'white'),
                edgecolor=style['marker_color'],
                linewidth=1.5,
                alpha=0.9,
                zorder=5
            )
            main_ax.add_patch(hex_patch)
        else:
            # Draw circle for other interactions
            marker_circle = Circle(
                (midpoint_x, midpoint_y), marker_radius,
                facecolor=style.get('marker_bg', 'white'),
                edgecolor=style['marker_color'],
                linewidth=1.5,
                alpha=0.9,
                zorder=5
            )
            main_ax.add_patch(marker_circle)
        
        # Add interaction symbol with dynamic font sizing
        font_size = 9 if text_length <= 1 else 8  # Smaller font for longer text
        
        text = main_ax.text(
            midpoint_x, midpoint_y,
            style['marker_text'],
            ha='center', va='center',
            fontsize=font_size, color=style['marker_color'],
            fontweight='bold',
            zorder=6
        )
        
        # Add subtle shadow/outline to the text for better readability
        if style.get('glow', False) and self.use_enhanced_styling:
            text.set_path_effects([
                path_effects.withStroke(linewidth=1.5, foreground='white'),
                path_effects.Normal()
            ])
    
    # Create legend for interaction types
    legend_title = "Interacting structural groups"
    legend_elements = []

    # Add residue example for legend
    residue_patch = Rectangle((0, 0), 1, 1, facecolor='white', 
                                edgecolor='black', label='Interacting structural groups')
    legend_elements.append(residue_patch)
    
    # Interaction type markers for legend
    for int_type, style in interaction_styles.items():
        # Only include interaction types that are present
        if self.interactions[int_type]:
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
    
    # Add solvent accessibility indicator
    if self.solvent_accessible:
        solvent_color = self.colors.get('solvent_color', '#ADD8E6')
        solvent_patch = Rectangle((0, 0), 1, 1, facecolor=solvent_color, 
                                    alpha=0.3, edgecolor=None, label='Solvent accessible')
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
