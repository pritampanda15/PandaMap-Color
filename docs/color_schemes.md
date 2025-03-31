# PandaMap Color Schemes Reference

PandaMap offers several predefined color schemes for visualizing protein-ligand interactions. This document provides details on each color scheme and how to create custom schemes.

## Built-in Color Schemes

PandaMap includes the following built-in color schemes:

### 1. Default

A balanced color scheme with good contrast and vibrant colors. Suitable for most visualizations.

- Background: Light blue-white (#F8F8FF)
- Ligand background: Light blue (#ADD8E6)
- Element colors:
  - Carbon: Dark grey (#444444)
  - Nitrogen: Rich blue (#3050F8)
  - Oxygen: Bright red (#FF2010)
  - Sulfur: Bright yellow (#FFFF30)
  - Phosphorus: Bright orange (#FF8000)
  - etc.

### 2. Colorblind

A carefully selected palette that remains distinguishable for people with various types of color vision deficiency.

- Background: White (#FFFFFF)
- Ligand background: Very light blue (#E6F5F9)
- Element colors:
  - Carbon: Gray (#666666)
  - Nitrogen: Blue (#0072B2)
  - Oxygen: Vermilion (#D55E00)
  - Sulfur: Yellow (#F0E442)
  - etc.

### 3. Monochrome

A grayscale palette for publications where color isn't available or desired.

- Background: White (#FFFFFF)
- Ligand background: Light gray (#F0F0F0)
- Element colors use various shades of gray to provide contrast

### 4. Dark

A dark-mode theme with light elements on a dark background. Excellent for presentations or reduced eye strain.

- Background: Almost black (#1A1A1A)
- Ligand background: Dark blue-gray (#2C3E50)
- Element colors use brighter versions to stand out against the dark background

### 5. Publication

A clean, professional style optimized for scientific publications. Uses subdued colors with good contrast.

- Background: White
- Ligand background: Very light gray (#F0F0F0)
- Element colors use medium-intensity tones that reproduce well in print

## Color Scheme Components

Each color scheme consists of the following components:

### Element Colors

Colors for different atom types in the ligand structure:

```json
"element_colors": {
    "C": "#444444",  // Carbon
    "N": "#3050F8",  // Nitrogen
    "O": "#FF2010",  // Oxygen
    "S": "#FFFF30",  // Sulfur
    "P": "#FF8000",  // Phosphorus
    "F": "#90E050",  // Fluorine
    "Cl": "#1FF01F", // Chlorine
    "Br": "#A62929", // Bromine
    "I": "#940094",  // Iodine
    "H": "#FFFFFF"   // Hydrogen
}
```

### Residue Colors

Colors for different types of amino acid residues:

```json
"residue_colors": {
    "hydrophobic": "#FFD700",    // Gold for hydrophobic
    "polar": "#00BFFF",          // Deep sky blue for polar
    "charged_pos": "#FF6347",    // Tomato for positive
    "charged_neg": "#32CD32",    // Lime green for negative
    "aromatic": "#DA70D6",       // Orchid for aromatic
    "default": "white"           // Default white
}
```

### Interaction Styles

Styles for different types of protein-ligand interactions:

```json
"interaction_styles": {
    "hydrogen_bonds": {
        "color": "#2E8B57",       // Line color
        "linestyle": "-",         // Solid line
        "linewidth": 1.8,         // Line thickness
        "marker_text": "H",       // Text in the marker
        "marker_color": "#2E8B57", // Color of the text
        "marker_bg": "#E0FFE0",   // Background color of marker
        "name": "Hydrogen Bond",  // Name for legend
        "glow": true              // Whether to add glow effect
    },
    // Other interaction types follow the same pattern
}
```

### Background Colors

Overall background styling:

```json
"background_color": "#F8F8FF",  // Main background
"ligand_bg_color": "#ADD8E6",   // Ligand area background
"solvent_color": "#ADD8E6"      // Solvent accessibility indicator
```

## Creating Custom Color Schemes

You can create custom color schemes in two ways:

### 1. Python Dictionary

Define your color scheme as a Python dictionary:

```python
from pandamap import ProtLigMapper

my_scheme = {
    "element_colors": {
        "C": "#333333",
        "N": "#5555FF",
        "O": "#FF5555",
        # Add other elements as needed
    },
    "residue_colors": {
        "hydrophobic": "#FFAA00",
        "polar": "#55AAFF",
        # Add other residue types as needed
    },
    "interaction_styles": {
        "hydrogen_bonds": {
            "color": "#00AA00",
            "linestyle": "-",
            "linewidth": 2.0,
            "marker_text": "H",
            "marker_color": "#00AA00",
            "marker_bg": "#CCFFCC",
            "name": "Hydrogen Bond",
            "glow": True
        },
        # Add other interaction types as needed
    },
    "background_color": "#FFFFFF",
    "ligand_bg_color": "#EEEEFF",
    "solvent_color": "#EEEEFF"
}

# Use your custom scheme
mapper = ProtLigMapper(
    "complex.pdb",
    color_scheme=my_scheme
)
```

### 2. JSON File

Create a JSON file with your color scheme:

```json
{
    "element_colors": {
        "C": "#333333",
        "N": "#5555FF",
        "O": "#FF5555"
    },
    "residue_colors": {
        "hydrophobic": "#FFAA00",
        "polar": "#55AAFF"
    },
    "interaction_styles": {
        "hydrogen_bonds": {
            "color": "#00AA00",
            "linestyle": "-",
            "linewidth": 2.0,
            "marker_text": "H",
            "marker_color": "#00AA00",
            "marker_bg": "#CCFFCC",
            "name": "Hydrogen Bond",
            "glow": true
        }
    },
    "background_color": "#FFFFFF",
    "ligand_bg_color": "#EEEEFF",
    "solvent_color": "#EEEEFF"
}
```

Then load it:

```python
mapper = ProtLigMapper(
    "complex.pdb",
    color_scheme="path/to/my_scheme.json"
)
```

## Tips for Creating Effective Color Schemes

1. **Maintain good contrast** between interaction lines and the background
2. **Use visually distinct colors** for different interaction types
3. **Consider color blindness** - test your scheme with tools like [Color Oracle](https://colororacle.org/)
4. **Keep print reproduction in mind** if designing for publications
5. **Be consistent with scientific conventions** where they exist (e.g., oxygen in red)
