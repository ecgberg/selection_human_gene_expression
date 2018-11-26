source activate python3.4
snakemake --rerun-incomplete --keep-going -j 500 --latency-wait 120 --cluster "sbatch --account pritch --job-name {params.job_name} --mem {params.memory} --time {params.run_time} -o {params.error_out_file}.out -e {params.error_out_file}.error --cpus-per-task {params.cores}"
