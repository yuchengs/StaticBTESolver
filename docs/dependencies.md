# Dependencies

StaticBTESolver can be compiled to a CPU version, or a GPU-accelerated version if you
have CUDA support. To use the CPU version, the only requirement is `Petsc`. To use the
GPU version, you need at least one CUDA-enabled GPU card and `CUDAToolkit` installed.
Detailed dependencies for both versions are summarized below.

## CPU version

### `PETSc`

Please consult PETSc documentation [here](https://www.mcs.anl.gov/petsc/documentation/installation.html).

### `FLTK` (optional)

gmsh's GUI uses FLTK. So if you want to visualize the mesh generated, you have to install FLTK first.

```
  curl -O https://www.fltk.org/pub/fltk/1.3.5/fltk-1.3.5-source.tar.gz
  tar zxvf fltk-1.3.5-source.tar.gz
  cd fltk-1.3.5
  ./configure
  make && sudo make install
```
  
### `OpenCASCADE` (optional)

OpenCASCADE is a free CAD kernel. If gmsh's build-in kernel is enough for you, you may skip this dependency. However,
we may support `STEP` file format in the future, which will probably use `OpenCASCADE`.

### `gmsh` (optional)

gmsh is a free mesh generator. To install, consult gmsh documentation [here](https://gmsh.info). Since we use gmsh api for C++ internally, you will need to 
compile gmsh as a shared library.

```
  git clone http://gitlab.onelab.info/gmsh/gmsh.git
  cd gmsh
  mkdir build
  cmake -S. -Bbuild -DENABLE_BUILD_DYNAMIC=1
  cd build
  make && sudo make install
```

## GPU version

### `CUDAToolkit` (required)

Please consult documentation [here](https://developer.nvidia.com/cuda-toolkit). We recommend
you to use the latest version.

### `viennacl` (required)

ViennaCL is a GPU-accelerated linear algebra library. It is a header library, so you
do not need to worry about how to install it.

### `openmpi` (required)

You may also use `mpich`.