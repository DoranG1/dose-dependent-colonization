#!/bin/bash
#SBATCH --job-name=submitPermuteCategoryPhylogeny
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --mail-user=dorang@stanford.edu
#SBATCH -p relman

source activate ~/miniconda3/envs/snakemakeR4.0

Rscript permuteCategoryPhylogenyOC.R
