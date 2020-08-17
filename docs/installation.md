# Installation

The library use C++17 features, so you should use a compiler supporting that.
Also, the minimum required version for `cmake` is 3.17.

Clone this repository:
```
git clone --recurse-submodules https://github.com/yuchengs/StaticBTESolver.git
cd StaticBTESolver
```
Then try to generate building tools with
```
mkdir build
cmake -S. -Bcpu-build -DCMAKE_BUILD_TYPE=Release # for CPU
cmake -S. -Bgpu-build -DCMAKE_BUILD_TYPE=Release -DUSE_GPU=on # for GPU
```
To use GPU version, you should have CUDA support. In case cmake cannot find path to nvcc,
provide the path for cmake:
```
cmake -S. -Bgpu-build -DCMAKE_BUILD_TYPE=Release -DUSE_GPU=on -DCMAKE_CUDA_COMPILER=path/to/nvcc
```
To build:
```
cd cpu-build # cd gpu-build for gpu version
make
```