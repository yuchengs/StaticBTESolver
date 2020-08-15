//
// Created by Yucheng Shi on 7/12/20.
//

#include "StaticBTESolver/BTEGeometry.h"
#include "gmsh.h"

BTEGeometry::BTEGeometry(const std::string& fileName, double L_x, double L_y, double L_z) {
    auto esc = fileName.find_last_of('.');
    std::string name = fileName.substr(0, esc);
    std::string ext = fileName.substr(esc + 1);
    assert(ext == "geo");

    this->L_x = L_x;
    this->L_y = L_y;
    this->L_z = L_z;

    mesh = new BTEMesh();
    mesh->L_x = L_x;
    mesh->L_y = L_y;
    mesh->L_z = L_z;

    gmsh::initialize();
//    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);
    int dim = gmsh::model::getDimension();
    mesh->dim = dim;
    gmsh::model::mesh::generate(dim);

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags,
                                coord,
                                parametricCoord,
                                -1,
                                -1,
                                false,
                                false);

    mesh->meshPts.reserve(nodeTags.size());
    for (int i = 0; i < nodeTags.size(); i++) {
        auto ptr = std::make_shared<staticbtesolver::Point>(coord[3 * i] * L_x, coord[3 * i + 1] * L_y, coord[3 * i + 2] * L_z);
        mesh->meshPts.push_back(ptr);
    }

    std::unordered_map<std::size_t, std::size_t> nodeTag_to_mesh_index;
    for (int i = 0; i < nodeTags.size(); i++) {
        nodeTag_to_mesh_index[nodeTags[i]] = i;
    }

    gmsh::vectorpair dimTags1D;
    gmsh::model::getEntities(dimTags1D, 1);

    for (auto p : dimTags1D) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags1D;
        gmsh::model::mesh::getElements(elementTypes,
                                       elementTags,
                                       nodeTags1D,
                                       p.first,
                                       p.second);
        for (int elementType_index = 0; elementType_index < elementTypes.size(); elementType_index++) {
            // dummy outer iteration: only allow segments
            assert(elementTypes[elementType_index] == 1);
            mesh->elements1D.reserve(elementTags[elementType_index].size());
            for (int element_index = 0; element_index < elementTags[elementType_index].size(); element_index++) {
                int start = nodeTag_to_mesh_index[nodeTags1D[elementType_index][2 * element_index]];
                int end = nodeTag_to_mesh_index[nodeTags1D[elementType_index][2 * element_index + 1]];
                auto ptr = std::make_shared<staticbtesolver::Segment>(start, end, p.second);
                mesh->elements1D.push_back(ptr);
            }
        }
    }

    gmsh::vectorpair dimTags2D;
    gmsh::model::getEntities(dimTags2D, 2);

    for (auto p : dimTags2D) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags2D;
        gmsh::model::mesh::getElements(elementTypes,
                                       elementTags,
                                       nodeTags2D,
                                       p.first,
                                       p.second);
        for (int elementType_index = 0; elementType_index < elementTypes.size(); elementType_index++) {
            // dummy outer iteration: only allow segments
            assert(elementTypes[elementType_index] == 2);
            mesh->elements2D.reserve(elementTags[elementType_index].size());
            for (int element_index = 0; element_index < elementTags[elementType_index].size(); element_index++) {
                int tri1 = nodeTag_to_mesh_index[nodeTags2D[elementType_index][3 * element_index]];
                int tri2 = nodeTag_to_mesh_index[nodeTags2D[elementType_index][3 * element_index + 1]];
                int tri3 = nodeTag_to_mesh_index[nodeTags2D[elementType_index][3 * element_index + 2]];
                auto ptr = std::make_shared<staticbtesolver::Triangle>(tri1, tri2, tri3, p.second);
                mesh->elements2D.push_back(ptr);
            }
        }
    }

    gmsh::vectorpair dimTags3D;
    gmsh::model::getEntities(dimTags3D, 3);

    for (auto p : dimTags3D) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags3D;
        gmsh::model::mesh::getElements(elementTypes,
                                       elementTags,
                                       nodeTags3D,
                                       p.first,
                                       p.second);
        for (int elementType_index = 0; elementType_index < elementTypes.size(); elementType_index++) {
            // dummy outer iteration: only allow segments
            assert(elementTypes[elementType_index] == 4);
            mesh->elements3D.reserve(elementTags[elementType_index].size());
            for (int element_index = 0; element_index < elementTags[elementType_index].size(); element_index++) {
                int tet1 = nodeTag_to_mesh_index[nodeTags3D[elementType_index][4 * element_index]];
                int tet2 = nodeTag_to_mesh_index[nodeTags3D[elementType_index][4 * element_index + 1]];
                int tet3 = nodeTag_to_mesh_index[nodeTags3D[elementType_index][4 * element_index + 2]];
                int tet4 = nodeTag_to_mesh_index[nodeTags3D[elementType_index][4 * element_index + 3]];
                auto ptr = std::make_shared<staticbtesolver::Tetrahedron>(tet1, tet2, tet3, tet4);
                mesh->elements3D.push_back(ptr);
            }
        }
    }
    if (dim == 2) {
//        std::cout << "Please specify boundary conditions for boundary # ";
        boundary_indices.reserve(dimTags1D.size());
        for (auto p : dimTags1D) {
//            std::cout << p.second << " ";
            boundary_indices.push_back(p.second);
        }
//        std::cout << std::endl;
    }
    if (dim == 3) {
//        std::cout << "Please specify boundary conditions for boundary # ";
        boundary_indices.reserve(dimTags2D.size());
        for (auto p : dimTags2D) {
//            std::cout << p.second << " ";
            boundary_indices.push_back(p.second);
        }
//        std::cout << std::endl;
    }
}

BTEMesh* BTEGeometry::export_mesh() {
    return mesh;
}

BTEGeometry::~BTEGeometry() {
    gmsh::finalize();
}
