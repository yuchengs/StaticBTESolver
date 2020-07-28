#!/usr/bin/env bash

# shellcheck disable=SC2046
cmake-build-debug/src/BTEcmd -x 1e-8 -y 1e-8 -z 1e-8 -t 4 -p 4 -w 1 -d 3 -I 1000 -g $(pwd)/tests/3DM3DG/mesh.mphtxt -b $(pwd)/tests/3DM3DG/inputbc.dat -m $(pwd)/tests/3DM3DG/band.dat > 3DM3DG.out
cmp --silent tests/3DM3DG/3DM3DG.sol 3DM3DG.out && echo '### 3DM3DG SUCCESS ###' || echo '### 3DM3DG ERROR ###'

# shellcheck disable=SC2046
cmake-build-debug/src/BTEcmd -x 1e-8 -y 1e-8 -z 0 -t 4 -p 4 -w 2 -d 3 -I 1000 -g $(pwd)/tests/3DM2DG/mesh.mphtxt -b $(pwd)/tests/3DM2DG/inputbc.dat -m $(pwd)/tests/3DM2DG/band.dat > 3DM2DG.out
cmp --silent tests/3DM2DG/3DM2DG.sol 3DM2DG.out && echo '### 3DM2DG SUCCESS ###' || echo '### 3DM2DG ERROR ###'

# shellcheck disable=SC2046
cmake-build-debug/src/BTEcmd -x 1e-8 -t 8 -p 1 -d 3 -I 1000 -g 100 -b $(pwd)/tests/3DM1DG/inputbc.dat -m $(pwd)/tests/3DM1DG/band.dat > 3DM1DG.out
cmp --silent tests/3DM1DG/3DM1DG.sol 3DM1DG.out && echo '### 3DM1DG SUCCESS ###' || echo '### 3DM1DG ERROR ###'

# shellcheck disable=SC2046
cmake-build-debug/src/BTEcmd -x 1e-8 -y 1e-8 -t 8 -p 8 -w 4pi -d 2 -I 10000 -g $(pwd)/tests/2DM2DG/mesh.mphtxt -b $(pwd)/tests/2DM2DG/inputbc.dat -m $(pwd)/tests/2DM2DG/band.dat > 2DM2DG.out
cmp --silent tests/2DM2DG/2DM2DG.sol 2DM2DG.out && echo '### 2DM2DG SUCCESS ###' || echo '### 2DM2DG ERROR ###'

# shellcheck disable=SC2046
cmake-build-debug/src/BTEcmd -x 1e-8 -t 8 -p 1 -d 2 -I 10000 -g 100 -b $(pwd)/tests/2DM1DG/inputbc.dat -m $(pwd)/tests/2DM1DG/band.dat > 2DM1DG.out
cmp --silent tests/2DM1DG/2DM1DG.sol 2DM1DG.out && echo '### 2DM1DG SUCCESS ###' || echo '### 2DM1DG ERROR ###'

rm *.out
