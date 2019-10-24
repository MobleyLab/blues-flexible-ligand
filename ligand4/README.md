# ligand 4
This directory provides the input files for running and analyzing simulations performed with ligand **4**. Note: for plotting the figures reported in the paper, first the scripts in this directory need to be run, followed by running the matlab scripts in the folder [`plots_paper`](../plots_paper).

## Contents

- [`MD-NCMC-flip`](MD-NCMC-flip): Directory containing input files for running MD/NCMCsimulations with 180 degree flip moves. The file [`run_MD_NCMC.sh`](MD-NCMC-flip/run_MD_NCMC.sh) contains commands used.
- [`MD`](MD): Directory containing input files for running 100 ns MD simulations. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`equi`](equi): Directory containing input files for equilibrating a solvated box of ligand **4** and c-Jun N-terminal kinase-1 inhibitor. The file [`run_equi.sh`](equi/run_equi.sh) contains commands used.
- complex_wat.prmtop - Amber compatible parameter file used to simulate ligand **4**
- ligand1.equi.rst - Coordinate file for an equilibrated solvated box of ligand **4** and c-Jun N-terminal kinase-1 inhibitor (PDB code: 2gmx)
