# Usage: bash run_NCMCstep_var.sh
# This scripts is for running MD/NCMC simulations with different total number of steps. Used for optimizing the NCMC protocol.

for noStep in 30 300 600 1000 1500 2050 3000
do

  # create folder and input scripts for each NCMC protocl
  mkdir move-${noStep}NCMC
  cd move-${noStep}NCMC

     # preparing and running simulations for each parameter
     sed "s/XYZ/$noStep/g" ../example.py > example.py
     sed "s/XYZ/$noStep/g" ../blues.pbs > blues.pbs 
 
     # run blues job
     python3 blues.py > blues.log
     
     # post-processing blues MD trajectory for torsion distribution
     cpptraj -i ../dihedral.ptraj
     
     # find accepted moves
     python3 ../../../scripts/findAcceptedMove.py gmx.log
     
     # calculate number of moves accepted as a function of total number of iterations proposed
     python3 ../../../scripts/findAcceptanceIteration.py gmx.log acc_ncmc_XYZNCMC.txt

  cd ..

done 
