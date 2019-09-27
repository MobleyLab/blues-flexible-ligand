# Usage: bash run_MD_NCMC.sh 

# running MD-NCMC via BLUES
python3 blues.py > blues.log

# get work distribution
grep log_ncmc gmx.log > work
cut -f 6 -d ' ' work > work_ncmc_noMove_alchRegion.txt
