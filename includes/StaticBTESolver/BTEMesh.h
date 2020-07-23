//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_BTEMESH_H
#define BTESOLVER_BTEMESH_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <cassert>
#include "utility.h"


class BTEMesh {
public:
    int dim{};
    double L_x, L_y, L_z;
    std::vector<std::shared_ptr<Point>> meshPts;
    std::vector<std::shared_ptr<Segment>> elements1D;
    std::vector<std::shared_ptr<Triangle>> elements2D;
    std::vector<std::shared_ptr<Tetrahedron>> elements3D;
    BTEMesh() = default;
    BTEMesh(std::ifstream& inFile, double L_x, double L_y = 0, double L_z = 0);
    BTEMesh(int N_cell, double L_x);
    BTEMesh(const BTEMesh& mesh) = delete;
};

#endif //BTESOLVER_BTEMESH_H
