#running MD-NCMC via BLUES
python3 blues.py > blues.log
# post-processing blues MD trajectory for torsion distribution
cpptraj -i dihedral.ptraj
