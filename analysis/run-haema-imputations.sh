#!/bin/bash
#SBATCH -J haema-imps
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=all
#SBATCH --time=0-23:00:00
#SBATCH --mem-per-cpu=10G


# Purge
module purge

# Load module
module load statistical/R/4.2.1

# Run imps
Rscript ./analysis/illustrative-example.R
