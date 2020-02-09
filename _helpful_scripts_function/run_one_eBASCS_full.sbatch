#!/bin/bash
#
#SBATCH -p stats          # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node
#SBATCH -t 0-4:00         # Running time
#SBATCH --mem 5000        # Memory request

module load R/3.4.2-fasrc01
module load gcc/7.1.0-fasrc01
module load python

export R_LIBS_USER=~/apps/R:~/apps/R_3.4.2_os7:$R_LIBS_USER


R CMD BATCH --no-save --no-restore "--args file.name='${1}'" run_eBASCS_full_fixedk.R ./sim_out/${1}_eBASCS.out 
