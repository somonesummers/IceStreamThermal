#!/bin/bash
#
#SBATCH --job-name=shearmarginStreamD
#SBATCH --output=matlab_test.“%j”.out
#SBATCH --error=matlab_test.“%j”.err
#SBATCH --partition=serc
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --mail-type=ALL
module load matlab
matlab -nodisplay < StreamDCopy.m
