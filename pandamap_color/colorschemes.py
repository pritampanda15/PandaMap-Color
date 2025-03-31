"""
Color schemes for PandaMap-Color.

This module contains pre-defined color schemes for protein-ligand interaction visualization.
"""

import json 

# Predefined color schemes
COLOR_SCHEMES = {
    'default': {
        'element_colors': {
            'C': '#444444',  # Dark Grey
            'N': '#3050F8',  # Rich blue
            'O': '#FF2010',  # Bright red
            'S': '#FFFF30',  # Bright yellow
            'P': '#FF8000',  # Bright orange
            'F': '#90E050',  # Bright green
            'Cl': '#1FF01F', # Green
            'Br': '#A62929', # Brown
            'I': '#940094',  # Purple
            'H': '#FFFFFF'   # White
        },
        'residue_colors': {
            'hydrophobic': '#FFD700',    # Gold for hydrophobic
            'polar': '#00BFFF',          # Deep sky blue for polar
            'charged_pos': '#FF6347',    # Tomato for positive
            'charged_neg': '#32CD32',    # Lime green for negative
            'aromatic': '#DA70D6',       # Orchid for aromatic
            'default': 'white'           # Default white
        },
        'interaction_styles': {
            'hydrogen_bonds': {
                'color': '#2E8B57',  # Sea Green
                'linestyle': '-',
                'linewidth': 1.8,
                'marker_text': 'H',
                'marker_color': '#2E8B57',
                'marker_bg': '#E0FFE0',  # Light green bg
                'name': 'Hydrogen Bond',
                'glow': True
            },
            'carbon_pi': {
                'color': '#4B4B4B',  # Darker gray
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'C-π',
                'marker_color': '#4B4B4B',
                'marker_bg': '#F8F8F8',  # Off-white
                'name': 'Carbon-Pi interaction',
                'glow': False
            },
            'pi_pi_stacking': {
                'color': '#8A2BE2',  # Blue Violet
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'π-π',
                'marker_color': '#8A2BE2',
                'marker_bg': '#F5E1FF',  # Very light purple
                'name': 'Pi-Pi stacking',
                'glow': True
            },
            'donor_pi': {
                'color': '#FF1493',  # Deep Pink
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'D',
                'marker_color': '#FF1493',
                'marker_bg': '#FFECF5',  # Very light pink
                'name': 'Donor-Pi interaction',
                'glow': False
            },
            'amide_pi': {
                'color': '#8B0000',  # Dark Red
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'A',
                'marker_color': '#8B0000',
                'marker_bg': '#FFE4E1',  # Misty Rose
                'name': 'Amide-Pi interaction',
                'glow': False
            },
            'hydrophobic': {
                'color': '#696969',  # Dim Gray
                'linestyle': ':',
                'linewidth': 1.5,
                'marker_text': 'h',
                'marker_color': '#696969',
                'marker_bg': '#F5F5F5',  # Very light gray
                'name': 'Hydrophobic',
                'glow': False
            }
        },
        'background_color': '#F8F8FF',  # Light background
        'ligand_bg_color': '#ADD8E6',   # Light blue for ligand
        'solvent_color': '#ADD8E6',     # Light blue for solvent accessibility
    },
    
    'colorblind': {
        'element_colors': {
            'C': '#666666',  # Gray
            'N': '#0072B2',  # Blue
            'O': '#D55E00',  # Vermilion
            'S': '#F0E442',  # Yellow
            'P': '#CC79A7',  # Pink
            'F': '#009E73',  # Green
            'Cl': '#009E73', # Green
            'Br': '#882255', # Purple-pink
            'I': '#882255',  # Purple-pink
            'H': '#FFFFFF'   # White
        },
        'residue_colors': {
            'hydrophobic': '#F0E442',    # Yellow
            'polar': '#56B4E9',          # Light blue
            'charged_pos': '#D55E00',    # Vermilion
            'charged_neg': '#009E73',    # Green
            'aromatic': '#CC79A7',       # Pink
            'default': 'white'           # Default white
        },
        'interaction_styles': {
            'hydrogen_bonds': {
                'color': '#009E73',  # Green
                'linestyle': '-',
                'linewidth': 2.0,
                'marker_text': 'H',
                'marker_color': '#009E73',
                'marker_bg': '#E6F5F0',
                'name': 'Hydrogen Bond',
                'glow': True
            },
            'carbon_pi': {
                'color': '#666666',  # Gray
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'C',
                'marker_color': '#666666',
                'marker_bg': '#F8F8F8',
                'name': 'Carbon-Pi interaction',
                'glow': False
            },
            'pi_pi_stacking': {
                'color': '#882255',  # Purple-pink
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'π',
                'marker_color': '#882255',
                'marker_bg': '#F5E6EB',
                'name': 'Pi-Pi stacking',
                'glow': True
            },
            'donor_pi': {
                'color': '#CC79A7',  # Pink
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'D',
                'marker_color': '#CC79A7',
                'marker_bg': '#F9EFF4',
                'name': 'Donor-Pi interaction',
                'glow': False
            },
            'amide_pi': {
                'color': '#D55E00',  # Vermilion
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'A',
                'marker_color': '#D55E00',
                'marker_bg': '#FAECE6',
                'name': 'Amide-Pi interaction',
                'glow': False
            },
            'hydrophobic': {
                'color': '#666666',  # Gray
                'linestyle': ':',
                'linewidth': 1.5,
                'marker_text': 'h',
                'marker_color': '#666666',
                'marker_bg': '#F5F5F5',
                'name': 'Hydrophobic',
                'glow': False
            }
        },
        'background_color': '#FFFFFF',  # White background
        'ligand_bg_color': '#E6F5F9',   # Very light blue for ligand
        'solvent_color': '#E6F5F9',     # Very light blue for solvent accessibility
    },
    
    'monochrome': {
        'element_colors': {
            'C': '#444444',  # Dark Gray
            'N': '#777777',  # Medium Gray
            'O': '#222222',  # Very Dark Gray
            'S': '#666666',  # Medium-Dark Gray
            'P': '#555555',  # Medium-Dark Gray
            'F': '#888888',  # Light Gray
            'Cl': '#888888', # Light Gray
            'Br': '#333333', # Dark Gray
            'I': '#333333',  # Dark Gray
            'H': '#AAAAAA'   # Very Light Gray
        },
        'residue_colors': {
            'hydrophobic': '#DDDDDD',    # Light Gray
            'polar': '#BBBBBB',          # Medium Light Gray
            'charged_pos': '#999999',    # Medium Gray
            'charged_neg': '#999999',    # Medium Gray
            'aromatic': '#BBBBBB',       # Medium Light Gray
            'default': 'white'           # White
        },
        'interaction_styles': {
            'hydrogen_bonds': {
                'color': '#444444',
                'linestyle': '-',
                'linewidth': 2.0,
                'marker_text': 'H',
                'marker_color': '#444444',
                'marker_bg': '#F5F5F5',
                'name': 'Hydrogen Bond',
                'glow': True
            },
            'carbon_pi': {
                'color': '#777777',
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'C',
                'marker_color': '#777777',
                'marker_bg': '#F5F5F5',
                'name': 'Carbon-Pi interaction',
                'glow': False
            },
            'pi_pi_stacking': {
                'color': '#555555',
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'π',
                'marker_color': '#555555',
                'marker_bg': '#F5F5F5',
                'name': 'Pi-Pi stacking',
                'glow': True
            },
            'donor_pi': {
                'color': '#666666',
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'D',
                'marker_color': '#666666',
                'marker_bg': '#F5F5F5',
                'name': 'Donor-Pi interaction',
                'glow': False
            },
            'amide_pi': {
                'color': '#888888',
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'A',
                'marker_color': '#888888',
                'marker_bg': '#F5F5F5',
                'name': 'Amide-Pi interaction',
                'glow': False
            },
            'hydrophobic': {
                'color': '#999999',
                'linestyle': ':',
                'linewidth': 1.5,
                'marker_text': 'h',
                'marker_color': '#999999',
                'marker_bg': '#F5F5F5',
                'name': 'Hydrophobic',
                'glow': False
            }
        },
        'background_color': '#FFFFFF',  # White background
        'ligand_bg_color': '#F0F0F0',   # Light gray for ligand
        'solvent_color': '#E6E6E6',     # Light gray for solvent accessibility
    },
    
    'dark': {
        'element_colors': {
            'C': '#E0E0E0',  # Light gray
            'N': '#7CB9E8',  # Brighter blue
            'O': '#FF6B6B',  # Brighter red
            'S': '#FFF07C',  # Brighter yellow
            'P': '#FFB347',  # Brighter orange
            'F': '#90EE90',  # Light green
            'Cl': '#90EE90', # Light green
            'Br': '#D8A0D8', # Light purple
            'I': '#B19CD8',  # Light indigo
            'H': '#FFFFFF'   # White
        },
        'residue_colors': {
            'hydrophobic': '#FFD700',    # Gold for hydrophobic
            'polar': '#00BFFF',          # Deep sky blue for polar
            'charged_pos': '#FF6347',    # Tomato for positive
            'charged_neg': '#7FFF00',    # Chartreuse for negative
            'aromatic': '#EE82EE',       # Violet for aromatic
            'default': '#E0E0E0'         # Light gray
        },
        'interaction_styles': {
            'hydrogen_bonds': {
                'color': '#50C878',  # Emerald green
                'linestyle': '-',
                'linewidth': 2.0,
                'marker_text': 'H',
                'marker_color': '#50C878',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Hydrogen Bond',
                'glow': True
            },
            'carbon_pi': {
                'color': '#D3D3D3',  # Light gray
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'C-π',
                'marker_color': '#D3D3D3',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Carbon-Pi interaction',
                'glow': False
            },
            'pi_pi_stacking': {
                'color': '#DA70D6',  # Orchid purple
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'π-π',
                'marker_color': '#DA70D6',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Pi-Pi stacking',
                'glow': True
            },
            'donor_pi': {
                'color': '#FFA07A',  # Light salmon
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'D',
                'marker_color': '#FFA07A',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Donor-Pi interaction',
                'glow': False
            },
            'amide_pi': {
                'color': '#FF6347',  # Tomato
                'linestyle': '--',
                'linewidth': 1.8,
                'marker_text': 'A',
                'marker_color': '#FF6347',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Amide-Pi interaction',
                'glow': False
            },
            'hydrophobic': {
                'color': '#A9A9A9',  # Dark gray
                'linestyle': ':',
                'linewidth': 1.5,
                'marker_text': 'h',
                'marker_color': '#A9A9A9',
                'marker_bg': '#1A1A1A',  # Very dark gray
                'name': 'Hydrophobic',
                'glow': False
            }
        },
        'background_color': '#1A1A1A',  # Almost black background
        'ligand_bg_color': '#2C3E50',   # Dark blue-gray for ligand
        'solvent_color': '#34495E',     # Slightly lighter blue-gray for solvent
    },
    
    'publication': {
        'element_colors': {
            'C': '#333333',  # Dark Gray
            'N': '#0066CC',  # Medium Blue
            'O': '#CC0000',  # Medium Red
            'S': '#CCCC00',  # Medium Yellow
            'P': '#CC6600',  # Medium Orange
            'F': '#00CC00',  # Medium Green
            'Cl': '#00CC00', # Medium Green
            'Br': '#996633', # Medium Brown
            'I': '#660066',  # Medium Purple
            'H': '#CCCCCC'   # Light Gray
        },
        'residue_colors': {
            'hydrophobic': '#E6C300',    # Gold
            'polar': '#0099CC',          # Blue
            'charged_pos': '#CC3300',    # Red
            'charged_neg': '#339900',    # Green
            'aromatic': '#993399',       # Purple
            'default': 'white'           # White
        },
        'interaction_styles': {
            'hydrogen_bonds': {
                'color': '#339900',
                'linestyle': '-',
                'linewidth': 1.5,
                'marker_text': 'H',
                'marker_color': '#339900',
                'marker_bg': 'white',
                'name': 'Hydrogen Bond',
                'glow': False
            },
            'carbon_pi': {
                'color': '#666666',
                'linestyle': '--',
                'linewidth': 1.2,
                'marker_text': 'C',
                'marker_color': '#666666',
                'marker_bg': 'white',
                'name': 'Carbon-Pi interaction',
                'glow': False
            },
            'pi_pi_stacking': {
                'color': '#993399',
                'linestyle': '--',
                'linewidth': 1.2,
                'marker_text': 'π',
                'marker_color': '#993399',
                'marker_bg': 'white',
                'name': 'Pi-Pi stacking',
                'glow': False
            },
            'donor_pi': {
                'color': '#CC3300',
                'linestyle': '--',
                'linewidth': 1.2,
                'marker_text': 'D',
                'marker_color': '#CC3300',
                'marker_bg': 'white',
                'name': 'Donor-Pi interaction',
                'glow': False
            },
            'amide_pi': {
                'color': '#990000',
                'linestyle': '--',
                'linewidth': 1.2,
                'marker_text': 'A',
                'marker_color': '#990000',
                'marker_bg': 'white',
                'name': 'Amide-Pi interaction',
                'glow': False
            },
            'hydrophobic': {
                'color': '#999999',
                'linestyle': ':',
                'linewidth': 1.0,
                'marker_text': 'h',
                'marker_color': '#999999',
                'marker_bg': 'white',
                'name': 'Hydrophobic',
                'glow': False
            }
        },
        'background_color': 'white',     # White background
        'ligand_bg_color': '#F0F0F0',    # Very light gray for ligand
        'solvent_color': '#F0F0F0',      # Very light gray for solvent
    }
}

def load_custom_color_scheme(json_file):
    """
    Load a custom color scheme from a JSON file.
    
    Parameters:
    -----------
    json_file : str
        Path to the JSON file with color scheme definition
        
    Returns:
    --------
    dict : Custom color scheme dictionary
    """
    try:
        with open(json_file, 'r') as f:
            custom_scheme = json.load(f)
        return custom_scheme
    except Exception as e:
        print(f"Error loading color scheme from {json_file}: {e}")
        print("Using default color scheme instead.")
        return COLOR_SCHEMES['default']
