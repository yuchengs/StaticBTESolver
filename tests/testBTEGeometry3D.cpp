//
// Created by Yucheng Shi on 7/12/20.
//

#include <vector>
#include "BTESolver/BTEGeometry.h"

using namespace std;

int main() {
    const string fileName = "cube.geo";
    BTEGeometry geo(fileName, 1, 1, 1);
    BTEMesh* mesh = geo.export_mesh();
}