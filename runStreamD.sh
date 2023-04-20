#!/bin/bash
#
#SBATCH --job-name=shearmarginStreamD
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=normal
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < StreamDCopy.m
