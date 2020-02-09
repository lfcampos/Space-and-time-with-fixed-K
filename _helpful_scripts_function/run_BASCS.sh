# loop over all files in the directory, print the filename
# (so we have some visual progress indicator), then submit the
# gzip jobs to SLURM
#
for FILE in ./sim_data/*.Rda; do
 echo ${FILE}
 sbatch -o ./sim_out/${FILE}_BASCS.stdout.txt -e ./sim_out/${FILE}_BASCS.stderr.txt run_one_BASCS.sbatch ${FILE}
 sleep 1 # pause to be kind to the scheduler
done

