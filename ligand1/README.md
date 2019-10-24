# ligand 1
This directory provides the input files for running and analyzing simulations performed with ligand **1**. Note: for plotting the figures reported in the paper, first the scripts in this directory need to be run, followed by running the matlab scripts in the folder [`plots_paper`](../plots_paper).

## Contents

- [`MD-NCMC-flip`](MD-NCMC-flip): Directory containing input files for running MD/NCMCsimulations with 180 degree flip moves. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`MD-NCMC-noMove`](MD-NCMC-noMove): Directory containing input files for running NCMC simulations with different alchemical regions. We do not propose any move and only alchemically switch the ligand off and on in the bidning pocket. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`MD-NCMC`](MD-NCMC): Directory containing input files for running MD/NCMC simulations with random torsional moves of selected rotatable bond. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`MD`](MD): Directory containing input files for running 100 ns MD simulations. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`NCMCstep_var`](NCMCstep_var): Directory containing input files for varying the number of NCMC steps to optimize NCMC protocol. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- [`equi`](equi): Directory containing input files for equilibrating a solvated box of ligand **1** and c-Jun N-terminal kinase-1 inhibitor. The file [`run_MD.sh`](MD/run_MD.sh) contains commands used.
- complex_wat.prmtop - Amber compatible parameter file used to simulate ligand **1**
- ligand1.equi.rst - Coordinate file for an equilibrated solvated box of ligand **1** and c-Jun N-terminal kinase-1 inhibitor (PDB code: 2gmx)


