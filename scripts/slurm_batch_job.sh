#!/bin/bash
#SBATCH -o /data/hpcdata/users/michai/slurm_job_output/job.%j.%N.out
#SBATCH -D /data/hpcdata/users/michai/
#SBATCH -J TEST                           # Specify a useful name for the job
#SBATCH --mem=2gb                             # Request 2GB of memory for the job
#SBATCH --time=00:02:00                       # Maximum job time of 12 hours
#SBATCH --cpus-per-task=16                    # Request 16 logical cpu cores

# Pick the partition and the matching billing account.
#SBATCH --partition=medium
#SBATCH --account=medium

# Enable software modules
source /etc/profile.d/modules.sh

# Load software modules you need
module add hpc/intel/2017
module add hpc/netcdf/intel/4.4.1.1

# Print date, run program, print date
date
hostname
date
