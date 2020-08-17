# Examples

## Examples using `BTEcmd`
   
### `2DM1DG`

Make sure you use files in `tests/2DM1DG`,
   
```
./BTEcmd -x 1e-8 -t 8 -p 1 \
   -d 2 -I 1000 \
   -g 100 \
   -b path/to/tests/2DM1DG/inputbc.dat \
   -m path/to/tests/2DM1DG/band.dat     
```
     
### `3DM1DG`

Make sure you use files in `tests/3DM1DG`,
   
```
./BTEcmd -x 1e-8 -t 8 -p 1 \
   -d 3 -I 1000 \
   -g 100 \
   -b path/to/tests/3DM1DG/inputbc.dat \
   -m path/to/tests/3DM1DG/band.dat     
```
     
### `2DM2DG`

Make sure you use files in `tests/2DM2DG`,
   
```
./BTEcmd -x 1e-4 -y 1e-4 -z 0 -d 2
   -t 8 -p 8 -I 1000000 
   -b path/to/tests/2DM2DG/inputbc.dat 
   -g path/to/tests/2DM2DG/mesh.mphtxt 
   -m path/to/tests/2DM2DG/band.dat
```
Note this configuration is slow, maybe try `x`, `y` with 1e-8 for faster convergence.
   
### `3DM2DG`

Make sure you use files in `tests/3DM2DG`,
   
```
./BTEcmd -x 1e-8 -y 1e-8 -z 0 -t 4 -p 4 \
   -d 3 -I 1000 \
   -g path/to/tests/3DM2DG/mesh.mphtxt \
   -b path/to/tests/3DM2DG/inputbc.dat \
   -m path/to/tests/3DM2DG/band.dat 
```
   
### `3DM3DG`

Make sure you use files in `tests/3DM3DG`,
   
```
./BTEcmd -x 1e-8 -y 1e-8 -z 1e-8 -t 4 -p 4 \
   -d 3 -I 1000 \
   -g path/to/tests/3DM3DG/mesh.mphtxt \
   -b path/to/tests/3DM3DG/inputbc.dat \
   -m path/to/tests/3DM3DG/band.dat 
```

### `2DM1DG`

Make sure you use files in `tests/2DM1DG`,

```
./BTEcmd -x 1e-8 -t 8 -p 1 \
    -d 2 -I 1000 \
    -g 100 \
    -b path/to/tests/2DM1DG/inputbc.dat \
    -m path/to/tests/2DM1DG/band.dat     
```

### `3DM1DG`

Make sure you use files in `tests/3DM1DG`,

```
./BTEcmd -x 1e-8 -t 8 -p 1 \
    -d 3 -I 1000 \
    -g 100 \
    -b path/to/tests/3DM1DG/inputbc.dat \
    -m path/to/tests/3DM1DG/band.dat     
```

### `2DM2DG`

Make sure you use files in `tests/2DM2DG`,

``` 
./BTEcmd -x 1e-4 -y 1e-4 -z 0 -d 2
    -t 8 -p 8 -I 1000000 
    -b path/to/tests/2DM2DG/inputbc.dat 
    -g path/to/tests/2DM2DG/mesh.mphtxt 
    -m path/to/tests/2DM2DG/band.dat
```
Note this configuration is slow, maybe try `x`, `y` with 1e-8 for faster convergence.

### `3DM2DG`

Make sure you use files in `tests/3DM2DG`,

```
./BTEcmd -x 1e-8 -y 1e-8 -z 0 -t 4 -p 4 \
    -d 3 -I 1000 \
    -g path/to/tests/3DM2DG/mesh.mphtxt \
    -b path/to/tests/3DM2DG/inputbc.dat \
    -m path/to/tests/3DM2DG/band.dat 
```
  
### `3DM3DG`

Make sure you use files in `tests/3DM3DG`,

```
./BTEcmd -x 1e-8 -y 1e-8 -z 1e-8 -t 4 -p 4 \
    -d 3 -I 1000 \
    -g path/to/tests/3DM3DG/mesh.mphtxt \
    -b path/to/tests/3DM3DG/inputbc.dat \
    -m path/to/tests/3DM3DG/band.dat 
```
  
## Examples using the API
  
``` cpp
#include <iostream>
#include <string>
#include "StaticBTESolver/StaticBTESolver.h"

using namespace std;

int main() {
    double L_x, L_y, L_z, Tref, errnum, qflux;
    int maxIter;
    bool backup;
    int ntheta, nphi;
    
    ifstream infile("inputdata3D.dat");
    string line;
    char newline;
    getline(infile, line);
    infile >> L_x >> newline;
    getline(infile, line);
    infile >> L_y >> newline;
    getline(infile, line);
    infile >> L_z >> newline;
    getline(infile, line);
    infile >> Tref >> newline;
    getline(infile, line);
    infile >> errnum >> newline;
    getline(infile, line);
    infile >> ntheta >> newline;
    getline(infile, line);
    infile >> nphi >> newline;
    getline(infile, line);
    infile >> maxIter >> newline;
    getline(infile, line);
    infile >> backup >> newline;
    infile.close();

    ifstream inFile1("Input-dispersion-relation-fp.dat");
    auto bands = new BTEBand(inFile1);
    ifstream inFile2("inputbc3D.dat");
    auto bcs = new BTEBoundaryCondition(inFile2);
    ifstream inFile3("mesh3D.mphtxt");
    auto mesh3D = new  BTEMesh(inFile3, L_x, L_y, L_z);

    StaticBTESolver solver(mesh3D, bcs, bands);
    solver.setParam(3, ntheta, nphi, Tref);
    solver.solve(maxIter);

    return 0;
}
```

## On SJTU PI cluster

The following script can be used to submit a job using two cores, each for
one GPU card.
``` bash
#!/bin/bash

#SBATCH --job-name=name
#SBATCH --partition=dgx2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=00-00:05:00
#SBATCH -n 2
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node 2
#SBATCH --mail-type=end
#SBATCH --mail-user=youremail

module purge
module load gcc/8.3.0-gcc-4.8.5
module load openmpi/3.1.5-gcc-8.3.0
module load cuda/10.1.243-gcc-8.3.0

srun --mpi=pmi2 gpu-build/src/BTEcmd -x 1e-8 -y 1e-8 \
           -t 8 -p 8 -d 2 -I 1000 \
           -g /path/to/2D_78000.mphtxt \
           -b /path/to/inputbc_2D.dat \
           -m /path/to/Input-dispersion-relation-fp.dat
```