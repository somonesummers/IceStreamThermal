#!/bin/bash
#
#SBATCH --job-name=shearmarginTest
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=normal
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < Site2.m
