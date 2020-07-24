//
// Created by Yucheng Shi on 7/12/20.
//

#ifndef BTESOLVER_BTEGEOMETRY_H
#define BTESOLVER_BTEGEOMETRY_H

#include <string>
#include <cassert>
#include <unordered_map>
#include "BTEMesh.h"

class BTEGeometry {
    double L_x, L_y, L_z;
    BTEMesh* mesh;
public:
    std::vector<int> boundary_indices;
    BTEGeometry() = delete;
    BTEGeometry(const std::string& fileName, double L_x, double L_y = 0, double L_z = 0);
    BTEMesh* export_mesh();
    ~BTEGeometry();
};

#endif //BTESOLVER_BTEGEOMETRY_H
