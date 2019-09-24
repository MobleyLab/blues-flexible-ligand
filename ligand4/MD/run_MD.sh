module load cuda/8.0.44
# running 100 ns MD on OpenMM
python3 mdin/md.py > PBS/md.log
# post-processing blues MD trajectory for torsion distribution
cpptraj -i dihedral.ptraj
