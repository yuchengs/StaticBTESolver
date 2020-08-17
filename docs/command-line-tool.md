# Command line tool

To save you from using the API (or any C++ programming), we provide a command line
utility `BTEcmd`.

You can find `BTEcmd` under `/cpu-build/src` or `/gpu-build/src`. This is the entry point to use
`StaticBTESolver` without having to compile and work with APIs.

## Arguments
Several arguments should be provided. Detailed input file format can be found in the 
API documentation.

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

## Sample Input Files

Some sample input files are in the `/tests` directory. For example, you may try

``` 
path/to/BTEcmd -g path/to/mesh2D.mphtxt \ 
    -m path/to/Input-dispersion-relation-fp.dat \ 
    -b path/to/Inputbc2D.dat \ 
    -d 2 -t 8 -p 8 \ 
    -x 1e-7 -y 1e-7 -z 1e-7
```

If you have installed the GPU version, use `mpirun`:

``` 
mpirun -np 2 path/to/BTEcmd -g path/to/mesh2D.mphtxt \
                 -m path/to/Input-dispersion-relation-fp.dat \
                 -b path/to/Inputbc2D.dat \
                 -d 2 -t 8 -p 8 \
                 -x 1e-7 -y 1e-7 -z 1e-7
```

