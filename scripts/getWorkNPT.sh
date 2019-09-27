for i in 10 100 500 #1000 2500 5000 10000 15000 
do
  echo $i
  cd noMove-${i}NCMC
  grep log_ncmc gmx.log > work
  cut -f 6 -d ' ' work > work_ncmc_noMove_${i}.txt
  mv work_ncmc_noMove_${i}.txt ../transfer/
  cd ../
done
