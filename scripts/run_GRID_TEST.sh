#!/bin/bash
##SBATCH -e slurm-%j.err
##SBATCH -o slurm-%j.out
##SBATCH --nodes=1

#SBATCH -o /data/hpcdata/users/michai/slurm_job_output/job.%j.%N.out

#SBATCH --cpus-per-task=16
#SBATCH --job-name=GRIDTEST

#SBATCH --time=00:02:00
#SBATCH --partition=short
#SBATCH --account=short
#SBATCH --chdir=/users/michai/mitgcmPy

#source /etc/profile.d/modules.sh

#module load python

python GRID_TEST.py

date
