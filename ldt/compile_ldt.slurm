#!/bin/bash -l
#SBATCH -t 00:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -A lp_ees_swm_ls_001
#SBATCH --cluster=genius
#SBATCH --mail-type FAIL,BEGIN,END
#SBATCH --mail-user michel.bechtold@kuleuven.be
#SBATCH -o compile_ldt_log.txt
#SBATCH -e compile_ldt_out.txt

####
# for larger lis simulations that need more than 6 cores you should add #SBATCH --exclude=r23i13n23

cd $SLURM_SUBMIT_DIR
source KUL_LDT_modules
./compile

