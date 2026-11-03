#!/bin/bash
#SBATCH --job-name=rclone
#SBATCH --account pedrini.edoardo
#SBATCH --mem=8GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=4  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="rcloneGdrivePersonal_data.err"
#SBATCH --output="rcloneGdrivePersonal_data.out"

echo "my job strart now" > rcloneGdrivePersonal_data.log;

date >> rcloneGdrivePersonal_data.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_rclone;

rclone copy -PL /beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/261030_BulkRNAseq_MG_Aletta/data \
gdrive_edo:work/analysis/HSR/absinta/bulkrnaseq/261030_BulkRNAseq_MG_Aletta/data

date >> rcloneGdrivePersonal_data.log;