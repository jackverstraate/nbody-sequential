#!/bin/bash
#SBATCH --job-name=nbody_bench
#SBATCH --partition=Centaurus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --output=slurm-%j.out

module purge
module load gcc

make clean
make

echo "Solar (N=10):"
(/usr/bin/time -p ./nbody solar.tsv 200 5000000 5000000 /dev/null) 2>&1

echo "N=100:"
(/usr/bin/time -p ./nbody 100 1 10000 10000 /dev/null) 2>&1

echo "N=1000 (calibration 200 steps):"
(/usr/bin/time -p ./nbody 1000 1 200 200 /dev/null) 2>&1