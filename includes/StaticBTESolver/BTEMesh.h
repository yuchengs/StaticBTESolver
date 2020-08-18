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
    double L_x{}, L_y{}, L_z{};
    std::vector<std::shared_ptr<staticbtesolver::Point>> meshPts;
    std::vector<std::shared_ptr<staticbtesolver::Segment>> elements1D;
    std::vector<std::shared_ptr<staticbtesolver::Triangle>> elements2D;
    std::vector<std::shared_ptr<staticbtesolver::Tetrahedron>> elements3D;
    BTEMesh() = default;

    /**
     * Constructor for 2D and 3D mesh
     * @param inFile    inFile should be an opened file stream of gmsh native format .msh
     * @param L_x       x-axis scale
     * @param L_y       y-axis scale
     * @param L_z       z-axis scale
     */
    BTEMesh(std::ifstream& inFile, double L_x, double L_y = 0, double L_z = 0);

    /**
     * Constructor for 1D mesh
     * @param N_cell
     * @param L_x
     */
    BTEMesh(int N_cell, double L_x);
    BTEMesh(const BTEMesh& mesh) = delete;
};

#endif //BTESOLVER_BTEMESH_H
