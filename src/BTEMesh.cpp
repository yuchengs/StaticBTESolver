//
// Created by Yucheng Shi on 7/6/20.
//

#include "StaticBTESolver/BTEMesh.h"

BTEMesh::BTEMesh(std::ifstream& inFile, double L_x, double L_y, double L_z) {
    if (!inFile.is_open()) {
        std::cout << "DEBUG: file not open" << std::endl;
        exit(1);
    }
    this->L_x = L_x;
    this->L_y = L_y;
    this->L_z = L_z;
    std::string line;
    while (getline(inFile, line)) {
        if (line == "4 Mesh # class") {
            break;
        }
    }
    getline(inFile, line);
    assert(line == "4 # version");

    getline(inFile, line);
    std::istringstream ss(line);
    ss >> this->dim;
    getline(inFile, line);
    ss.str(line);
    int num_meshPts;
    ss >> num_meshPts;
    meshPts.reserve(num_meshPts);

    while (getline(inFile, line)) {
        if (line == "# Mesh point coordinates") {
            break;
        }
    }

    for (int i = 0; i < num_meshPts; i++) {
        auto ptr = std::make_shared<Point>();
        if (this->dim == 3) {
            double x, y, z;
            inFile >> x >> y >> z;
            ptr->x = x * L_x;
            ptr->y = y * L_y;
            ptr->z = z * L_z;
        }
        else if (this->dim == 2) {
            double x, y;
            inFile >> x >> y;
            ptr->x = x * L_x;
            ptr->y = y * L_y;
        }
        else {
            inFile >> ptr->x;
            ptr->x *= L_x;
        }
        meshPts.push_back(std::move(ptr));
    }

    while (getline(inFile, line)) {
        if (line == "3 edg # type name") {
            break;
        }
    }
    // Populate 1D
    for (int i = 0; i < 4; i++)
        getline(inFile, line);
    ss.str(line);
    int num_1D;
    ss >> num_1D;
    elements1D.reserve(num_1D);
    getline(inFile, line);
    for (int i = 0; i < num_1D; i++) {
        auto ptr = std::make_shared<Segment>();
        inFile >> ptr->index[0] >> ptr->index[1];
        elements1D.push_back(std::move(ptr));
    }
    for (int i = 0; i < 4; i++)
        getline(inFile, line);
    for (int i = 0; i < num_1D; i++) {
        inFile >> elements1D[i]->entity_index;
    }

    if (this->dim == 1) return;
    while (getline(inFile, line)) {
        if (line == "3 tri # type name") {
            break;
        }
    }
    for (int i = 0; i < 4; i++) {
        getline(inFile, line);
    }
    ss.str(line);
    int num_2D;
    ss >> num_2D;
    elements2D.reserve(num_2D);
    getline(inFile, line);
    for (int i = 0; i < num_2D; i++) {
        auto ptr = std::make_shared<Triangle>();
        inFile >> ptr->index[0] >> ptr->index[1] >> ptr->index[2];
        elements2D.push_back(std::move(ptr));
    }

    if (this->dim == 2) return;
    for (int i = 0; i < 4; i++)
        getline(inFile, line);
    for (int i = 0; i < num_2D; i++) {
        inFile >> elements2D[i]->entity_index;
    }

    while (getline(inFile, line)) {
        if (line == "3 tet # type name") {
            break;
        }
    }
    for (int i = 0; i < 4; i++) {
        getline(inFile, line);
    }
    ss.str(line);
    int num_3D;
    ss >> num_3D;
    elements3D.reserve(num_3D);
    getline(inFile, line);
    for (int i = 0; i < num_3D; i++) {
        auto ptr = std::make_shared<Tetrahedron>();
        inFile >> ptr->index[0] >> ptr->index[1] >> ptr->index[2] >> ptr->index[3];
        elements3D.push_back(std::move(ptr));
    }
}

BTEMesh::BTEMesh(int N_cell, double L_x) {
    this->L_x = L_x;
    this->L_y = 0;
    this->L_z = 0;
    this->dim = 1;
    meshPts.reserve(N_cell + 1);
    for (int cell_index = 0; cell_index < N_cell + 1; cell_index++) {
        auto ptr = std::make_shared<Point>(L_x / N_cell * cell_index);
        meshPts.push_back(std::move(ptr));
    }
    elements1D.reserve(N_cell);
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        auto ptr = std::make_shared<Segment>(cell_index, cell_index + 1);
        elements1D.push_back(std::move(ptr));
    }
}



