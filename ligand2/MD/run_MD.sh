# Usage: bash run_MD.sh

module load cuda/8.0.44

mkdir rst PBS netcdf mdout mdinfo

# running 100 ns MD on OpenMM
python3 ../../scripts/mdin/md.py ../ligand2.equi.rst ../complex_wat.prmtop > PBS/md.log

# post-processing blues MD trajectory for torsion distribution
cpptraj -i dihedral.ptraj
