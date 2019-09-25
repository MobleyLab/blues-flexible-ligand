# blues-flexible-ligand
contains files for blues-flexible-ligand paper


# Input files for Stage 1a pose prediction subchallenge using HYBRID (OpenEye Scientific SOftware)

This directory (and its subdirectories) provides the input files of the receptors (reference pdb structures used for docking), the pdb codes of the respective reference ligands and the mol2 files of the docked poses of the 20 ligands that we prepared using HYBRID and submitted for Stage 1a pose prediction subchallenge. Here, we also provide all the scripts that we used to parametrize the receptors and the ligands using amber and tleap, respectively. Moreover, the scripts used to perform MD simulations are available in this directory.

## Contents

- [`ligand1`](ligand1): directory containing the mol2 files of the 2 selected docked poses of the ligand BACE_1 and the scripts used to run MD and MM-GBSA calculations on the rescpective protein-ligand complexes.
