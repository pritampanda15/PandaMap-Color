"""
Command-line interface for PandaMap-Color.
"""

import os
import sys
import argparse
from .pandamap import PandaMapColor
from .colorschemes import COLOR_SCHEMES, load_custom_color_scheme


def main():
    """
    Command-line interface for PandaMap-Color.
    """
    # Create an argument parser with a description
    parser = argparse.ArgumentParser(
        description="""
        PandaMap-Color: Visualize protein-ligand interactions from PDB files
        with customizable styling and coloring options.
        """)
    
    # Add arguments
    parser.add_argument('pdb_file', help='Path to PDB file with protein-ligand complex')
    
    parser.add_argument('--output', '-o', help='Output image file path')
    
    parser.add_argument('--ligand', '-l', help='Specific ligand residue name to analyze')
    
    parser.add_argument('--dpi', type=int, default=300, 
                      help='Image resolution (default: 300 dpi)')
    
    parser.add_argument('--title', '-t', help='Custom title for the visualization')
    
    parser.add_argument('--color-scheme', '-c', 
                      choices=list(COLOR_SCHEMES.keys()),
                      default='default',
                      help='Color scheme to use (default: default)')
    
    parser.add_argument('--custom-colors', 
                      help='Path to JSON file with custom color scheme')
    
    parser.add_argument('--simple-styling', action='store_true',
                      help='Use simple styling instead of enhanced effects')
    
    parser.add_argument('--no-color-by-type', action='store_true',
                      help='Disable coloring residues by amino acid type')
    
    parser.add_argument('--jitter', type=float, default=0.1,
                      help='Amount of positional randomization (0.0-1.0) for more natural look')
    
    parser.add_argument('--h-bond-cutoff', type=float, default=3.5,
                      help='Distance cutoff for hydrogen bonds in Angstroms')
    
    parser.add_argument('--pi-stack-cutoff', type=float, default=5.5,
                      help='Distance cutoff for pi-stacking interactions in Angstroms')
    
    parser.add_argument('--hydrophobic-cutoff', type=float, default=4.0,
                      help='Distance cutoff for hydrophobic interactions in Angstroms')
    
    parser.add_argument('--figsize', type=float, nargs=2, default=(12, 12),
                      help='Figure size in inches (width height)')
    
    # New arguments for enhanced functionality
    parser.add_argument('--generate-report', action='store_true',
                      help='Generate a detailed interaction report')
    
    parser.add_argument('--report-file',
                      help='Path for the interaction report (default: <pdb_name>_interactions_report.txt)')
    
    parser.add_argument('--use-dssp', action='store_true',
                      help='Use DSSP for solvent accessibility calculation if available')
    
    parser.add_argument('--show-directionality', action='store_true',
                      help='Show interaction directionality with arrows')
    
    parser.add_argument('--ionic-cutoff', type=float, default=4.0,
                      help='Distance cutoff for ionic interactions in Angstroms')
    
    parser.add_argument('--halogen-cutoff', type=float, default=3.5,
                      help='Distance cutoff for halogen bonds in Angstroms')
    
    # Parse arguments
    args = parser.parse_args()
    
    try:
        print("╔════════════════════════════════════════════════════╗")
        print("║                 PandaMap-Color                     ║")
        print("╚════════════════════════════════════════════════════╝")
        print(f"Processing: {os.path.basename(args.pdb_file)}")
        print(f"Color scheme: {args.color_scheme}")
        print(f"Style: {'Simple' if args.simple_styling else 'Enhanced'}")
        
        # Determine color scheme
        if args.custom_colors:
            color_scheme = load_custom_color_scheme(args.custom_colors)
            print(f"Using custom color scheme from: {args.custom_colors}")
        else:
            color_scheme = args.color_scheme
        
        # Initialize mapper
        mapper = PandaMapColor(
            args.pdb_file, 
            ligand_resname=args.ligand,
            color_scheme=color_scheme,
            use_enhanced_styling=not args.simple_styling
        )
        
        # Check if we should use the enhanced run_analysis method
        if hasattr(mapper, 'run_analysis') and callable(getattr(mapper, 'run_analysis')):
            # Use the enhanced run_analysis method with report generation
            print("Using enhanced analysis pipeline...")
            result = mapper.run_analysis(
                output_file=args.output,
                use_dssp=args.use_dssp,
                generate_report=args.generate_report,
                report_file=args.report_file
            )
        else:
            # Fall back to manual step-by-step approach
            print("Using standard analysis pipeline...")
            # Detect interactions with custom cutoffs
            print("Detecting interactions...")
            mapper.detect_interactions(
                h_bond_cutoff=args.h_bond_cutoff,
                pi_stack_cutoff=args.pi_stack_cutoff,
                hydrophobic_cutoff=args.hydrophobic_cutoff,
                ionic_cutoff=getattr(args, 'ionic_cutoff', 4.0),
                halogen_bond_cutoff=getattr(args, 'halogen_cutoff', 3.5)
            )
            
            # Estimate solvent accessibility
            print("Estimating solvent accessibility...")
            if args.use_dssp and hasattr(mapper, 'calculate_dssp_solvent_accessibility'):
                try:
                    mapper.calculate_dssp_solvent_accessibility()
                except Exception as e:
                    print(f"DSSP calculation failed: {e}. Falling back to estimation.")
                    mapper.estimate_solvent_accessibility()
            else:
                mapper.estimate_solvent_accessibility()
            
            # Generate visualization
            print("Generating visualization...")
            from .visualization import visualize
            result = visualize(
                mapper,
                output_file=args.output,
                figsize=args.figsize,
                dpi=args.dpi,
                title=args.title,
                color_by_type=not args.no_color_by_type,
                jitter=args.jitter,
                show_directionality=args.show_directionality if hasattr(args, 'show_directionality') else False
            )
            
            # Generate report if requested
            if args.generate_report and hasattr(mapper, 'generate_interaction_report'):
                report_file = args.report_file
                if report_file is None:
                    base_name = os.path.splitext(os.path.basename(args.pdb_file))[0]
                    report_file = f"{base_name}_interactions_report.txt"
                
                print(f"Generating interaction report: {report_file}")
                try:
                    mapper.generate_interaction_report(output_file=report_file)
                    print(f"Report saved to: {report_file}")
                except Exception as e:
                    print(f"Error generating report: {e}")
        
        # Check if visualization was successful
        if result and not (isinstance(result, str) and result.startswith("Error")):
            print("✓ Analysis complete!")
            print(f"✓ Visualization saved to: {result}")
            return 0
        else:
            print(f"✗ Visualization failed: {result}")
            return 1
        
    except Exception as e:
        print(f"✗ Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())