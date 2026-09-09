CHARMM-GUI Membrane Setup Extractor
A robust, standard-library Python command-line tool to recursively scan CHARMM-GUI project directories, auto-detect metadata, and generate a comprehensive structural and setup report for your molecular dynamics simulations.

Overview
When building membrane or protein-membrane systems using CHARMM-GUI, the critical design parameters (such as planned lipid ratios, intended water box thickness, orientation constraints, and salt concentrations) are scattered across multiple stream (.str), parameter (.prm), and coordinate (.gro) files.

This script searches your project directory, parses these files, cross-references planned builds with the final structural coordinate file (step5_input.gro), and compiles a clean, publication-ready text report.

Features
The script extracts, calculates, and reports:
System Dimensions: Specified and calculated box geometry, lattice vectors, angles (γ), area (nm^2 & "\unicode" 〖x212B〗^2), and volume in both Nanometers and Angstroms.
Lipid Composition & Leaflet Distribution: Planned vs. final counts for every lipid species in both the upper and lower leaflets, flagging any differences (deltas).
Component Classification: Categorized counts for lipids, water molecules, and ions parsed from the definitive final coordinate file.
Water Information: Water thickness/padding values alongside exact solvent counts.
Ion details: Declared types, counts, concentration parameters, and final confirmed ion counts.
Peptide/Protein PPM Data: Orientation, tilt angles, insertion/depth details, and cross-sectional or lipid-exclusion area coordinates.
Atom and Residue Inventories: Detailed breakdowns of every single residue type present in the GROMACS coordinate file.
Workflow & Command Extraction: MDP execution sequences and command lines found in build logs (README, shell scripts).
Target Files Auto-Detected
The script recursively scans your target folder to find:

input.config.dat (System build configurations)
step3_nlipids_upper.prm (Upper leaflet lipid plans)
step3_nlipids_lower.prm (Lower leaflet lipid plans)
step4_components.str (Component lists)
step4.3_ion.str (Ionization setups)
step2_orient.out (PPM peptide orientation outputs)
step5_input.gro (Final assembled GROMACS coordinate file)
Requirements
Python 3.9 or higher
No external dependencies (uses the Python Standard Library only)
Works on Windows, macOS, and Linux.
Usage
1. Basic Execution
Place the script in your main directory containing your CHARMM-GUI outputs (or point to it) and run:

Bash

python3 extract_membrane_setup_1.py
This scans the current directory (.) and outputs a file named system_information.txt.

2. Scanning a Specific Folder
If your simulation files are in a different directory:

Bash

python3 extract_membrane_setup_1.py /path/to/charmm-gui-outputs/
3. Output Customization
Change the name/destination of the generated report:

Bash

python3 extract_membrane_setup_1.py -o my_system_report.txt
4. Handling Custom Residues
If your system contains custom lipids, ions, or water models that are not automatically recognized by standard dictionaries, you can pass them as comma-separated lists:

Bash

python3 extract_membrane_setup_1.py \
  --lipid-resnames "MYLIP1,MYLIP2" \
  --ion-resnames "ZN,ZN2+" \
  --water-resnames "TIP4"
5. Disable GROMACS Run Template Generation
To skip generating the automated GROMACS run template sequence, use:

Bash

python3 extract_membrane_setup_1.py --no-run-template
CLI Options Help
To see all available options directly in your terminal, run:

Bash

python3 extract_membrane_setup_1.py --help
text

positional arguments:
  root                  Project directory to search recursively. (Default: ".")

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output report file path. (Default: "system_information.txt")
  --water-resnames WATER_RESNAMES
                        Additional water residue names (comma-separated).
  --ion-resnames ION_RESNAMES
                        Additional ion residue names (comma-separated).
  --lipid-resnames LIPID_RESNAMES
                        Additional lipid residue names (comma-separated).
  --no-run-template     Do not generate GROMACS execution template.
Typical Directory Layout
The tool handles nested folders beautifully. It will successfully find files in standard CHARMM-GUI layouts:

text

charmm-gui-project/
├── input.config.dat
├── step2_orient.out
├── step3_nlipids_lower.prm
├── step3_nlipids_upper.prm
├── step4.3_ion.str
├── step4_components.str
└── gromacs/
    ├── step5_input.gro
    ├── step6.1_equilibration.mdp
    ├── step6.2_equilibration.mdp
    └── ...
Understanding the Output Report
The output text file is divided into clear sections with visual borders:

Detected Input Files: Paths to all auto-located raw files.
Box and Solvation Settings: Water padding thickness, specified box parameters, and final box coordinates converted into both Nanometers and Angstroms.
Lipid Counts by Leaflet and Species: A tabular breakdown of planned vs. actual lipids per species per leaflet.
Confirmed Final Component Counts: Quantified values of lipids, waters, ions, and proteins.
Water Information: Confirmed water count and residue mappings.
Ion Information: Declared concentrations vs. actual ionic counts.
Peptide Orientation / PPM Information: Orientation parameters if a peptide/protein is present.
Peptide Exclusion Area: Cross-sectional footprint of the protein in the membrane.
Total Atom Count & Residue Breakdowns: Detailed index table of all residues.
Exact Final Box Vectors: Raw system matrix vectors.
Run Instructions / Step Sequence: Step-by-step shell commands auto-generated based on the sequence of your GROMACS .mdp files.
Warnings / Extraction Notes: Informational notes highlighting discrepancies (such as lipid count deviations or missing files).
