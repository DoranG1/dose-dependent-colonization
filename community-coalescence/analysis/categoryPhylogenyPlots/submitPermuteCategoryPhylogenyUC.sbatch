#!/bin/bash
#SBATCH --job-name=submitPermuteCategoryPhylogeny
#SBATCH --time=48:00:00
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=dorang@stanford.edu
#SBATCH -p relman

source activate ~/miniconda3/envs/snakemakeR4.0

Rscript permuteCategoryPhylogenyUC.R
