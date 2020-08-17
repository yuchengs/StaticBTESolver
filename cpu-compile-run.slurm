#!/bin/bash

#SBATCH --job-name=cpu-2dm2dg-case5
#SBATCH --partition=small
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=00-00:30:00
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --mail-type=end
#SBATCH --mail-user=youremail

module purge
module load gcc/8.3.0-gcc-4.8.5
module load openmpi/3.1.5-gcc-8.3.0

rm -rf cpu-build
cmake3 -S. -Bcpu-build -DCMAKE_BUILD_TYPE=Release
cd cpu-build
make
cd ..
cpu-build/src/BTEcmd -x 1e-8 -y 1e-8 \
           -t 8 -p 8 -d 2 -I 1000 \
           -g /lustre/home/acct-umjbh/umjbh-b/yuchengs/StaticBTESolver/tests/testcase/2Dcase5/2D_78000.mphtxt \
           -b /lustre/home/acct-umjbh/umjbh-b/yuchengs/StaticBTESolver/tests/testcase/2Dcase5/inputbc_2D.dat \
           -m /lustre/home/acct-umjbh/umjbh-b/yuchengs/StaticBTESolver/tests/testcase/2Dcase5/Input-dispersion-relation-fp.dat
