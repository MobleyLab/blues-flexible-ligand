# Usage: bash run_umbrella_sampling.sh 
# This script is for running and post-processing MD simulations for each window in umbrella sampling simulations.
# 152 windows, same force constant, variable force constant

mkdir PBS dcd csv rst dihed

# here I am doing 10 ns MD simulations for each window sequentially, can be easily parallelized and is recommended for faster turnaround time
umbNo=0
while read line;
do

  splitted=($line)
  arr=${splitted[0]}   # centers of each window
  K=${splitted[1]}     # force constant
  echo $umbNo
  # running 10ns MD simulation
  python3 ../simulate_umbrella.py $umbNo $arr $K > PBS/umbrella_${umbNo}.log  
  # calculating dihedral distribution
  sed "s/XYZ/XXX/g" ../dihedral.ptraj > dihedralXXX.ptraj
  cpptraj -i dihedralXXX.ptraj
  rm dihedralXXX.ptraj

  ((umbNo++))

done < centers.dat

# calculating PMF from the dihedral angle distribution of each window
python3 analyze_umbrella.py 152 > PMF_umbrella.txt
# deleting first line for matlab processing
sed -i '1d' PMF_umbrella.txt
 
