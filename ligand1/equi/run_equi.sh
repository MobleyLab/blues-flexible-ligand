# Usage: bash run_equi.sh

module load cuda/8.0.44

# runnning antechamber to get AM1-BCC partial charges for the ligand
antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -s 2  > antechamber.log

# parameterizing the ligand
parmchk -i LIG.mol2 -f mol2 -o LIG.frcmod
tleap -f ../../scripts/leap_LIG.in > leap_LIG.log

# creating simulation box with tleap
tleap -f leap_2gmx.in > leap_2gmx.log

# making directories for writing MD output
mkdir rst PBS netcdf mdout mdinfo

# running equilibration on OpenMM
for i in $(seq 1 8)
do
   echo $i
   python3 ../../scripts/mdin/step${i}.py > PBS/step${i}.log
   date
done

#copying final structure for MD and MD/NCMC simulation
cp rst/step8.rst.250000 ../ligand1.equi.rst
