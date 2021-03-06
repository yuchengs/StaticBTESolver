# StaticBTESolver

## Dependencies

1. `PETSc`: Please consult PETSc documentation [here](https://www.mcs.anl.gov/petsc/documentation/installation.html). Note we recommend 
you to download and install MPI and/or LAPACK (by providing command line arguments to `./configure`) even if you have already installed them on your machine. 
Instructions for your reference:
  ```
  ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
  make all test
  ```
  After installation, take a note on the `PETSC_DIR` and the `PETSC_ARCH` variables, you will need them later.
  
2. `FLTK` (optional): gmsh's GUI uses FLTK. So if you want to visualize the mesh generated, you have to install FLTK first.
  ```
  curl -O https://www.fltk.org/pub/fltk/1.3.5/fltk-1.3.5-source.tar.gz
  tar zxvf fltk-1.3.5-source.tar.gz
  cd fltk-1.3.5
  ./configure
  make && sudo make install
  ```
  
3. `OpenCASCADE` (optional): OpenCASCADE is a free CAD kernel. If gmsh's build-in kernel is enough for you, you may skip this dependency. However,
we may support `STEP` file format in the future, which will probably use `OpenCASCADE`.

4. `gmsh` (optional): gmsh is a free mesh generator. To install, consult gmsh documentation [here](gmsh.info). Since we use gmsh api for C++ internally, you will need to 
compile gmsh as a shared library.
  ```
  git clone http://gitlab.onelab.info/gmsh/gmsh.git
  cd gmsh
  mkdir build
  cmake -S. -Bbuild -DENABLE_BUILD_DYNAMIC=1
  cd build
  make && sudo make install
  ```
  
## Installation
```
git clone --recurse-submodules https://github.com/yuchengs/StaticBTESolver.git
cd StaticBTESolver
```
Find the `CMakeLists.txt` and set `PETSC_DIR` and `PETSC_ARCH` to correct ones. Then try to generate building tools with
```
mkdir build
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DPETSC_DIR=... -DPETSC_ARCH=... # for CPU
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DUSE_GPU=on # for GPU
```
To use GPU version, you should have CUDA support. In case cmake cannot find path to nvcc,
provide the path for cmake:
```
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DUSE_GPU=on -DCMAKE_CUDA_COMPILER=path/to/nvcc
```
To build:
```
cd build
make
```

For CPU version, you can try ctest in `build`.

## Running
After building, you can find a command line utility called `BTEcmd` under `/build/src`. This is the entry point to use
`StaticBTESolver` without having to compile and work with APIs.

Several arguments should be provided:

- `-g [string]` (required): provide a file with either gmsh native file format (with extension `.geo`) or 
    COMSOL native file format (with extension `.mphtxt`), or an integer specifying the number of cells if in the setting
    of 1DG. Note that a `.geo` file only specifies
    geometry information while a `.mphtxt` specifies mesh information. Internally, a `.geo` file is passed to 
    `BTEGeometry`, which meshes the geometry and exports to a `BTEMesh` object. In contrast, a `.mphtxt`
    file is passed to the constructor of `BTEMesh` directly. 
- `-m [string]` (required): provide a file specifying the bands information in the similar format of `Input-dispersion-relation-fp.dat`.
- `-b [string]` (depends): provide a file specifying the boundary conditions of each
    boundary. The file format is similar to `inputbc2D.dat`. If your geometry information is provided with a `.geo` file, this argument
    is optional. If absent, you will have to specify boundary conditions interactively.
    However, if you use a `.mphtxt` file, this argument is mandatory.
- `-d [integer]` (required): provide an integer specifying DM, __not__ DG. (DG is handled automatically)  
- `-t [integer]` (required): CADOM parameter. Number of discrete polar angle.
- `-p [integer]` (required): CADOM parameter. Number of discrete azimuthal angle.
- `-T [float]` (optional): reference temperature. Default at 300.
- `-I [integer]` (optional): maximum number of iteration. Default at 10000.
- `-x [float]` (optional): x-axis scale. Default at 0.
- `-y [float]` (optional): y-axis scale. Default at 0.
- `-z [float]` (optional): z-axis scale. Default at 0.

Some sample input files are in the `/tests` directory. For example, you may try
```$xslt
BTEcmd -g path/to/mesh2D.mphtxt \
    -m path/to/Input-dispersion-relation-fp.dat \
    -b path/to/Inputbc2D.dat \
    -d 2 -t 8 -p 8 \
    -x 1e-7 -y 1e-7 -z 1e-7
```
If you install the GPU version, use `mpirun`:
```
# 2 processes for 2 GPU
mpirun -np 2 ./BTEcmd -g path/to/mesh2D.mphtxt \
                 -m path/to/Input-dispersion-relation-fp.dat \
                 -b path/to/Inputbc2D.dat \
                 -d 2 -t 8 -p 8 \
                 -x 1e-7 -y 1e-7 -z 1e-7
```

## Example Usage and Testing

- `2DM1DG`: make sure you use files in `tests/2DM1DG`,
    ```$xslt
    ./BTEcmd -x 1e-8 -t 8 -p 1 \
        -d 2 -I 1000 \
        -g 100 \
        -b path/to/tests/2DM1DG/inputbc.dat \
        -m path/to/tests/2DM1DG/band.dat     
    ```
- `3DM1DG`: make sure you use files in `tests/3DM1DG`,
    ```$xslt
    ./BTEcmd -x 1e-8 -t 8 -p 1 \
        -d 3 -I 1000 \
        -g 100 \
        -b path/to/tests/3DM1DG/inputbc.dat \
        -m path/to/tests/3DM1DG/band.dat     
    ```
- `2DM2DG`: make sure you use files in `tests/2DM2DG`,
    ```
    ./BTEcmd -x 1e-4 -y 1e-4 -z 0 -d 2
        -t 8 -p 8 -I 1000000 
        -b path/to/tests/2DM2DG/inputbc.dat 
        -g path/to/tests/2DM2DG/mesh.mphtxt 
        -m path/to/tests/2DM2DG/band.dat
    ```
    Note this configuration is slow, maybe try `x`, `y` with 1e-8 for faster convergence.
- `3DM2DG`: make sure you use files in `tests/3DM2DG`,
    ```$xslt
    ./BTEcmd -x 1e-8 -y 1e-8 -z 0 -t 4 -p 4 \
        -d 3 -I 1000 \
        -g path/to/tests/3DM2DG/mesh.mphtxt \
        -b path/to/tests/3DM2DG/inputbc.dat \
        -m path/to/tests/3DM2DG/band.dat 
    ```
- `3DM3DG`: make sure you use files in `tests/3DM3DG`,
    ```$xslt
    ./BTEcmd -x 1e-8 -y 1e-8 -z 1e-8 -t 4 -p 4 \
        -d 3 -I 1000 \
        -g path/to/tests/3DM3DG/mesh.mphtxt \
        -b path/to/tests/3DM3DG/inputbc.dat \
        -m path/to/tests/3DM3DG/band.dat 
    ```
  
(__For Developers__): For all those examples, correct results are cached. So for example, you may compare your result from the correct one 
using the following for `3DM3DG`:
```$xslt
    ./BTEcmd -x 1e-8 -y 1e-8 -z 1e-8 -t 4 -p 4 \
        -d 3 -I 1000 \
        -g path/to/tests/3DM3DG/mesh.mphtxt \
        -b path/to/tests/3DM3DG/inputbc.dat \
        -m path/to/tests/3DM3DG/band.dat > outfile
    diff outfile path/to/tests/3DM3DG/3DM3DG.sol
```
If there is no output, you are fine. Otherwise your output is incorrect. On macOS, 
you can do batch checking by calling the shell script `test.sh` from the source directory.
On Linux, results may be different from macOS, so if you do not pass all tests from
linux, that does not mean there is any problem.

       
