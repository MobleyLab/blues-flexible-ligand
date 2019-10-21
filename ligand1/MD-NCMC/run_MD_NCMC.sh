# Usage: bash run_MD_NCMC.sh 

# running MD-NCMC via BLUES
python3 blues.py > blues.log

# post-processing blues MD trajectory for torsion distribution
cpptraj -i dihedral.ptraj

# find accepted moves
python3 ../../scripts/findAcceptedMove.py gmx.log 
