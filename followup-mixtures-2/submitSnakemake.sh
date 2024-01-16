#!/bin/bash
#SBATCH --job-name=submitSnakemake
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=kxue@stanford.edu
#SBATCH -p relman
source activate /home/users/kxue/relman/users/kxue/miniconda3/envs/snakemake-new
#snakemake -np --use-conda --cores 24 --cluster 'sbatch -t 500 --mem=96g -c 24 -p relman' -j 25 --max-jobs-per-second 3 --max-status-checks-per-second 3
snakemake --use-conda --cores 24 --cluster 'sbatch -t 96:00:00 --mem=96g -c 24 -p relman --output=slurm/%j.out' -j 25 --max-jobs-per-second 3 --max-status-checks-per-second 3
#snakemake -R mergeSNPsAllHouseholds --use-conda --cores 48 --cluster 'sbatch -t 1500 --mem=256g -c 48 -p relman'
