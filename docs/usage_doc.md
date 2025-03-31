# PandaMap-Color Usage Guide

This document provides detailed information on how to use PandaMap-Color for protein-ligand interaction visualization.

## Command Line Interface

PandaMap-Color provides a command-line interface for quick visualization:

```bash
pandamap-color [PDB_FILE] [OPTIONS]
```

### Required Arguments

- `pdb_file`: Path to PDB file with protein-ligand complex

### Options

- `--output`, `-o`: Output image file path (default: [pdb_name]_interactions.png)
- `--ligand`, `-l`: Specific ligand residue name to analyze
- `--dpi`: Image resolution (default: 300 dpi)
- `--title`, `-t`: Custom title for the visualization
- `--color-scheme`, `-c`: Color scheme to use (choices: default, colorblind, monochrome, dark, publication)
- `--custom-colors`: Path to JSON file with custom color scheme
- `--simple-styling`: Use simple styling instead of enhanced effects
- `--no-color-by-type`: Disable coloring residues by amino acid type
- `--jitter`: Amount of positional randomization (0.0-1.0) for more natural look
- `--h-bond-cutoff`: Distance cutoff for hydrogen bonds in Angstroms (default: 3.5)
- `--pi-stack-cutoff`: Distance cutoff for pi-stacking interactions in Angstroms (default: 5.5)
- `--hydrophobic-cutoff`: Distance cutoff for hydrophobic interactions in Angstroms (default: 4.0)
- `--figsize`: Figure size in inches (width height) (default: 12 12)

### Examples

```bash
# Basic usage
pandamap-color complex.pdb -o interactions.png

# Specify a ligand and title
pandamap-color complex.pdb -l LIG -t "Drug X with Target Y" -o drug_x_interactions.png

# Use colorblind-friendly scheme with higher resolution
pandamap-color complex.pdb -c colorblind --dpi 600 -o colorblind_interactions.png

# Custom hydrogen bond detection threshold
pandamap-color complex.pdb --h-bond-cutoff 3.2 -o custom_hbond_interactions.png
```

## Python API

For more control and integration with your workflow, use the Python API:

```python
from pandamap_color import PandaMapColor

# Initialize mapper
mapper = PandaMapColor(
    pdb_file="complex.pdb",
    ligand_resname="LIG",  # Optional: specify ligand residue name
    color_scheme="default",
    use_enhanced_styling=True
)

# Detect interactions with custom parameters
mapper.detect_interactions(
    h_bond_cutoff=3.5,
    pi_stack_cutoff=5.5,
    hydrophobic_cutoff=4.0
)

# Estimate solvent accessibility
mapper.estimate_solvent_accessibility()

# Generate visualization with custom parameters
mapper.visualize(
    output_file="interactions.png",
    figsize=(12, 12),
    dpi=300,
    title="My Protein-Ligand Complex",
    color_by_type=True,
    jitter=0.1
)
```

## Main Class: PandaMapColor

The `PandaMapColor` class is the main class for analyzing protein-ligand interactions.

### Constructor Parameters

- `pdb_file` (str): Path to the PDB file with protein and ligand
- `ligand_resname` (str, optional): Specific residue name of the ligand to focus on
- `color_scheme` (str or dict): Color scheme to use, either a key from COLOR_SCHEMES or a custom dictionary
- `use_enhanced_styling` (bool): Whether to use enhanced styling effects (gradients, shadows, etc.)

### Methods

- `detect_interactions(h_bond_cutoff=3.5, pi_stack_cutoff=5.5, hydrophobic_cutoff=4.0)`: Detect all interactions between protein and ligand
- `estimate_solvent_accessibility()`: Estimate which residues might be solvent accessible
- `visualize(output_file=None, figsize=(12, 12), dpi=300, title=None, color_by_type=True, jitter=0.1)`: Generate a 2D visualization of protein-ligand interactions
- `run_analysis(output_file=None)`: Run the complete analysis pipeline

## Customizing Visualizations

### Color Schemes

PandaMap-Color provides several built-in color schemes:

```python
from pandamap_color import PandaMapColor, COLOR_SCHEMES

# List available color schemes
print(COLOR_SCHEMES.keys())
# Output: dict_keys(['default', 'colorblind', 'monochrome', 'dark', 'publication'])

# Use a specific color scheme
mapper = PandaMapColor("complex.pdb", color_scheme="dark")
```

### Custom Color Schemes

You can define a custom color scheme as a Python dictionary or a JSON file:

```python
# Python dictionary
custom_colors = {
    "element_colors": {
        "C": "#444444",
        "N": "#3050F8",
        # more elements...
    },
    "residue_colors": {
        "hydrophobic": "#FFD700",
        "polar": "#00BFFF",
        # more residue types...
    },
    "interaction_styles": {
        "hydrogen_bonds": {
            "color": "#2E8B57",
            "linestyle": "-",
            "linewidth": 1.8,
            "marker_text": "H",
            "marker_color": "#2E8B57",
            "marker_bg": "#E0FFE0",
            "name": "Hydrogen Bond",
            "glow": True
        },
        # more interaction types...
    },
    "background_color": "#F8F8FF",
    "ligand_bg_color": "#ADD8E6",
    "solvent_color": "#ADD8E6"
}

# Use the custom scheme
mapper = PandaMapColor(
    pdb_file="complex.pdb",
    color_scheme=custom_colors
)
```

### JSON Color Scheme

You can also save your color scheme as a JSON file:

```python
import json

# Save color scheme to file
with open("my_scheme.json", "w") as f:
    json.dump(custom_colors, f, indent=2)

# Load and use the color scheme
mapper = PandaMapColor(
    pdb_file="complex.pdb",
    color_scheme="my_scheme.json"
)
```

## Advanced Usage

### Accessing Interaction Data

You can access the detected interactions directly:

```python
mapper = PandaMapColor("complex.pdb")
mapper.detect_interactions()

# Access hydrogen bonds
hbonds = mapper.interactions["hydrogen_bonds"]
for hbond in hbonds:
    lig_atom = hbond["ligand_atom"]
    prot_atom = hbond["protein_atom"]
    residue = hbond["protein_residue"]
    distance = hbond["distance"]
    
    print(f"H-bond: Ligand atom {lig_atom.get_name()} - Protein residue {residue.resname} {residue.id[1]}")
    print(f"         Protein atom: {prot_atom.get_name()}, Distance: {distance:.2f} Å")
```

### Handling Different File Formats

PandaMap-Color includes utilities for handling different file formats:

```python
from pandamap_color.utils import MultiFormatParser

# Create a parser that can handle multiple formats
parser = MultiFormatParser()

# Parse a structure (supports PDB, mmCIF, PDBQT)
structure = parser.parse_structure("complex.pdbqt")
```

### Integration with Other Tools

PandaMap-Color integrates well with other Python tools:

```python
import matplotlib.pyplot as plt
from pandamap_color import PandaMapColor

# Create a figure
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Create two mappers with different color schemes
mapper1 = PandaMapColor("complex.pdb", color_scheme="default")
mapper2 = PandaMapColor("complex.pdb", color_scheme="dark")

# Run analyses
mapper1.detect_interactions()
mapper1.estimate_solvent_accessibility()
mapper2.detect_interactions()
mapper2.estimate_solvent_accessibility()

# Generate visualizations on separate axes
mapper1.visualize(axes[0], title="Default Color Scheme")
mapper2.visualize(axes[1], title="Dark Color Scheme")

# Save the combined figure
plt.tight_layout()
plt.savefig("comparison.png", dpi=300)
```

## Troubleshooting

### Missing Ligand

If PandaMap-Color can't find your ligand:

```python
try:
    mapper = PandaMapColor("complex.pdb")
except ValueError as e:
    if "No ligand" in str(e):
        print("Try specifying the ligand name:")
        print("mapper = PandaMapColor('complex.pdb', ligand_resname='LIG')")
    else:
        raise e
```

### Visualization Quality

If the visualization appears cluttered:

```python
mapper = PandaMapColor("complex.pdb")
mapper.detect_interactions()
mapper.estimate_solvent_accessibility()

# Use a larger figure size and higher DPI
mapper.visualize(
    output_file="high_quality.png",
    figsize=(16, 16),  # Larger figure
    dpi=600,           # Higher resolution
    jitter=0.15        # More space between elements
)
```