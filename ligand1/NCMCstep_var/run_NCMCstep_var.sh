# Usage: bash run_NCMCstep_var.sh
# This scripts is for running MD/NCMC simulations with different total number of steps. Used for optimizing the NCMC protocol.

for noStep in 30 300 600 1000 1500 2050 3000
do

  # create folder and input scripts for each NCMC protocl
  mkdir move-${noStep}NCMC
  cd move-${noStep}NCMC
  sed "s/XYZ/$noStep/g" ../example.py > example.py
  sed "s/XYZ/$noStep/g" ../blues.pbs > blues.pbs 
 
  # run blues job
  qsub blues.pbs
  cd ..

done 
