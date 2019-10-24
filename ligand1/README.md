# ligand 1
This directory provides the input files for running and analyzing simulations performed with ligand **1**.

## Contents

- [`MD-NCMC-flip`](MD-NCMC-flip): directory containing input files for running MD/NCMCsimulations with 180 degree flip moves.
- [`MD-NCMC-noMove`](MD-NCMC-noMove): directory containing input files for running NCMC simulations with different alchemical regions. We do not propose any move and only alchemically switch the ligand off and on in the bidning pocket.
- [`MD-NCMC`](MD-NCMC): directory containing input files for running MD/NCMC simulations with random torsional moves of selected rotatable bond.
- [`MD`](MD): Directory containing input files for running 100 ns MD simulations. The file [`MD/run_MD.sh`](run_MD.sh)run_MD.sh contains commands used.
- [`NCMCstep_var`](NCMCstep_var): directory containing input files for varying the number of NCMC steps to optimize NCMC protocol.
- [`equi`](equi): directory containing input files for equilibrating a solvated box of ligand **1** and c-Jun N-terminal kinase-1 inhibitor.
- complex_wat.prmtop - Amber compatible parameter file used to simulate ligand **1**
- ligand1.equi.rst - Coordinate file for an equilibrated solvated box of ligand **1** and c-Jun N-terminal kinase-1 inhibitor (PDB code: 2gmx)


