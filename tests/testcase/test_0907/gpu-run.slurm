#!/bin/bash

#SBATCH --job-name=test
#SBATCH --partition=dgx2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=00-00:05:00
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node 1


module purge
module load gcc/8.3.0-gcc-4.8.5
module load openmpi/3.1.5-gcc-8.3.0
module load cuda/10.1.243-gcc-8.3.0

srun --mpi=pmi2 /lustre/home/acct-umjbh/umjbh-b/yuchengs/StaticBTESolver/gpu-build/src/BTEcmd -x 1e-9 -y 1e-9 \
           -t 1 -p 1 -d 2 -I 1000 \
           -g ./cc0.mphtxt \
           -b ./inputbc.dat \
           -m ./Input-dispersion-relation-fp.dat
