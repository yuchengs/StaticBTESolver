//
// Created by Yucheng Shi on 7/8/20.
//

#include "StaticBTESolver/StaticBTESolver.h"

StaticBTESolver::StaticBTESolver(BTEMesh* mesh, BTEBoundaryCondition* bcs, BTEBand* bands) {
    this->mesh = mesh;
    this->bcs = bcs;
    this->bands = bands;
}

void StaticBTESolver::setParam(int DM, int num_theta, int num_phi, double WFACTOR, double T_ref) {
    this->num_theta = num_theta;
    this->num_phi = num_phi;
    this->WFACTOR = WFACTOR;
    this->T_ref = T_ref;
    this->DM = DM;
    if (mesh->dim == 2) {
        N_cell = mesh->elements2D.size();
        if (DM == 3) {
            N_dir = 4 * num_theta * num_phi;
        }
        else if (DM == 2) {
            N_dir = 4 * num_phi;
        }
        N_face = 3;
    }
    else if (mesh->dim == 3) {
        N_cell = mesh->elements3D.size();
        if (DM == 3) {
            N_dir = 8 * num_theta * num_phi;
        }
        else if (DM == 2) {
            N_dir = 4 * num_theta * num_phi;
        }
        N_face = 4;
    }
    N_band = bands->size();
    cell_band_temperature = std::vector<std::vector<double>>(N_cell, std::vector<double>(N_band));
    cell_band_density = std::vector<std::vector<double>>(N_cell, std::vector<double>(N_band));
    ee_curr = vector3D<double>(N_cell, vector2D<double>(N_dir, std::vector<double>(N_band)));
    cell_temperature.resize(N_cell);
}

void StaticBTESolver::_get_cell_temperature(int band_index) {
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        double e0 = 0;
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
            e0 += ee_curr[cell_index][dir_index][band_index] * control_angles[dir_index];
        }
        cell_band_temperature[cell_index][band_index] = e0 / (*bands)[band_index].Ctot + T_ref;
    }
}

void StaticBTESolver::_recover_temperature() {
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        double RRn = 0, Rnn = 0;
        for (int band_index = 0; band_index < N_band; band_index++) {
            RRn += (*bands)[band_index].Lr * cell_band_temperature[cell_index][band_index];
            Rnn += (*bands)[band_index].Lr;
        }
        cell_temperature[cell_index] = RRn / Rnn;
    }
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        for (int band_index = 0; band_index < N_band; band_index++) {
            cell_band_density[cell_index][band_index] = (*bands)[band_index].Ctot / (4 * PI) * (cell_temperature[cell_index] - T_ref);
        }
    }
}

void StaticBTESolver::_get_const_coefficient() {
    dv_dot_normal_cache.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));
    S_dot_normal_cache.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));
    a_f_total.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));
    for (int band_index = 0; band_index < N_band; band_index++) {
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
            vector2D<double> Ke(N_cell, std::vector<double>(N_cell, 0));
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                for (int face_index = 0; face_index < N_face; face_index++) {
                    double area = cell_face_area[cell_index][face_index];
                    dv_dot_normal_cache[band_index][dir_index][cell_index].push_back(dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]));
                    S_dot_normal_cache[band_index][dir_index][cell_index].push_back(dot_prod(S[dir_index], cell_face_normal[cell_index][face_index]));
                    double temp = (*bands)[band_index].group_velocity * (*bands)[band_index].relaxation_time
                                  * S_dot_normal_cache[band_index][dir_index][cell_index][face_index] * area;
                    if (dv_dot_normal_cache[band_index][dir_index][cell_index][face_index] >= 0) {
                        Ke[cell_index][cell_index] += temp;
                    } else if (cell_neighbor_indices[cell_index][face_index] >= 0) {
                        Ke[cell_index][cell_neighbor_indices[cell_index][face_index]] += temp;
                    }
                    a_f_total[band_index][dir_index][cell_index].push_back(temp);
                }
                Ke[cell_index][cell_index] += cell_volume[cell_index] * control_angles[dir_index];
            }
            for (int i = 0; i < N_cell; i++) {
                for (int j = 0; j < N_cell; j++) {
                    if (Ke[i][j] != 0) {
                        Ke_serialized.push_back(band_index);
                        Ke_serialized.push_back(dir_index);
                        Ke_serialized.push_back(i);
                        Ke_serialized.push_back(j);
                        Ke_serialized.push_back(Ke[i][j]);
                    }
                }
            }
        }
    }
}

std::vector<double> StaticBTESolver::_get_coefficient(int dir_index, int band_index) {
    std::vector<double> Re(N_cell, 0);
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        for (int face_index = 0; face_index < N_face; face_index++) {
            for (auto bc : *bcs) {
                if (dv_dot_normal_cache[band_index][dir_index][cell_index][face_index] >= 0) continue;
                if (cell_neighbor_indices[cell_index][face_index] == bc.index) {
                    if (bc.type == 1) {
                        Re[cell_index] -= (*bands)[band_index].Ctot / (4 * PI)
                                          * (bc.temperature - T_ref)
                                          * a_f_total[band_index][dir_index][cell_index][face_index];
                    }
                    else if (bc.type == 2) {
                        double einsum = 0;
                        double temp = 0;
                        // CAUTION
                        for (int dir_index_inner = 0; dir_index_inner < N_dir; dir_index_inner++) {
                            if (dv_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index] > 0) {
                                if (mesh->dim == 3) {
                                    einsum += ee_prev[cell_index][dir_index_inner][band_index]
                                              * S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index]
                                              * control_angles[dir_index_inner];
                                    temp += S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index]
                                            * control_angles[dir_index_inner];
                                }
                                else if (mesh->dim == 2) {
                                    einsum += ee_prev[cell_index][dir_index_inner][band_index]
                                              * S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index];
                                }
                            }
                        }
                        if (mesh->dim == 3) einsum = einsum / temp;
                        else if (mesh->dim == 2) einsum = einsum / PI;
                        Re[cell_index] -= einsum * a_f_total[band_index][dir_index][cell_index][face_index];
                    }
                    else if (bc.type == 31) {
                        int itheta = dir_index / (4 * num_phi);
                        int iphi = dir_index % (4 * num_phi);
                        if (mesh->dim == 2) {
                            if (iphi >= 2 * num_phi) {
                                iphi = 6 * num_phi - iphi - 1;
                            } else {
                                iphi = 2 * num_phi - iphi - 1;
                            }
                        }
                        else if (mesh->dim == 3) {
                            iphi = 4 * num_phi - iphi - 1;
                        }
                        Re[cell_index] -= ee_prev[cell_index][itheta * 4 * num_phi + iphi][band_index]
                                          * a_f_total[band_index][dir_index][cell_index][face_index];
                    }
                    else if (bc.type == 32) {
                        int itheta = dir_index / (4 * num_phi);
                        int iphi = dir_index % (4 * num_phi);
                        if (mesh->dim == 2) {
                            iphi = 4 * num_phi - iphi - 1;
                            Re[cell_index] -= ee_prev[cell_index][itheta * 4 * num_phi + iphi][band_index]
                                              * a_f_total[band_index][dir_index][cell_index][face_index];
                        }
                        else if (mesh->dim == 3) {
                            int ix = iphi / (2 * num_phi);
                            int isx = iphi % (2 * num_phi);
                            isx = 2 * num_phi - isx - 1;
                            Re[cell_index] -= ee_prev[cell_index][itheta * 4 * num_phi + isx + ix * 2 * num_phi][band_index]
                                              * a_f_total[band_index][dir_index][cell_index][face_index];
                        }
                    }
                    else if (bc.type == 33) {
                        if (mesh->dim == 2) {
                            std::cout << "Boundary Condition Type " << bc.type << " not supported. Abort." << std::endl;
                            exit(1);
                        }
                        int iz = dir_index / (4 * num_phi);
                        int isz = dir_index % (4 * num_phi);
                        iz = 2 * num_phi - iz - 1;
                        Re[cell_index] -= ee_prev[cell_index][isz + iz * 4 * num_phi][band_index]
                                          * a_f_total[band_index][dir_index][cell_index][face_index];
                    }
                    else {
                        std::cout << "Boundary Condition Type " << bc.type << " not supported. Abort." << std::endl;
                        exit(1);
                    }
                }
            }
        }
        Re[cell_index] += cell_band_density[cell_index][band_index] * cell_volume[cell_index] * control_angles[dir_index];
    }
    return Re;
}

std::vector<double> StaticBTESolver::_solve_matrix(vector2D<double>& Ke, std::vector<double>& Re) {
    int* nnz = new int[N_cell];
    for (int i = 0; i < N_cell; i++) {
        *(nnz + i) = 0;
        for (int j = 0; j < N_cell; j++)
            if (Ke[i][j] != 0) {
                *(nnz + i) = *(nnz + i) + 1;
            }
    }

    static char help[] = "Solving matrix.\n\n";
    const double offset = 1e16;
    PetscInitialize(nullptr, nullptr, (char*)nullptr, help);
    Mat A;
    Vec xxx, bbb;
    PetscScalar v, u, yyy[N_cell];
    PetscInt iii, jjj, one, ix[N_cell];
    KSP ksp;
    MatCreateSeqAIJ(PETSC_COMM_SELF, N_cell, N_cell, 3, nnz, &A);
    for (iii = 0;iii < N_cell; iii++){
        for (jjj = 0; jjj < N_cell; jjj++){
            if (Ke[iii][jjj] != 0) {
                v = Ke[iii][jjj] * offset;
                MatSetValues(A, 1, &iii, 1, &jjj, &v, ADD_VALUES);
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecCreate(PETSC_COMM_WORLD, &bbb);
    VecSetSizes(bbb, PETSC_DECIDE, N_cell);
    VecSetFromOptions(bbb);
    VecDuplicate(bbb, &xxx);
    for (iii = 0; iii < N_cell; iii++){
        u = Re[iii] * offset;
        VecSetValues(bbb, 1, &iii, &u, ADD_VALUES);
    }

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);

    KSPSolve(ksp, bbb, xxx);

    one = N_cell;
    std::vector<double> x(N_cell);

    for (iii = 0; iii < N_cell; iii++)
        ix[iii] = iii;

    VecGetValues(xxx, one, ix, yyy);
    for (iii = 0; iii < N_cell; iii++)
        x[iii]=yyy[iii];
    delete[] nnz;
    VecDestroy(&xxx);
    VecDestroy(&bbb);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    return x;
}

double StaticBTESolver::_get_margin() {
    double margin = 0;
    for (int band_index = 0; band_index < N_band; band_index++){
        for (int dir_index = 0; dir_index < N_dir; dir_index++){
            for (int cell_index = 0; cell_index < N_cell; cell_index++)
                margin += pow((ee_curr[cell_index][dir_index][band_index] - ee_prev[cell_index][dir_index][band_index]), 2);
        }
    }
    return (margin / N_band / N_cell);
}

void StaticBTESolver::_get_heat_flux() {
    bc_band_heat_flux = vector2D<double>(bcs->size(), std::vector<double>(N_band, 0));
    bc_heat_flux = std::vector<double>(bcs->size(), 0);
    for (int bc_index = 0; bc_index < bcs->size(); bc_index++) {
        BoundaryCondition bc = (*bcs)[bc_index];
        for (int band_index = 0; band_index < N_band; band_index++) {
            double cell_length = 0;
            double heat = 0;
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                for (int face_index = 0; face_index < N_face; face_index++) {
                    if (cell_neighbor_indices[cell_index][face_index] == bc.index) {
                        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                            // TODO: check if control angle is needed
                            if (DM == 3 && mesh->dim == 3) {
                                heat += ee_curr[cell_index][dir_index][band_index]
                                        * dv_dot_normal_cache[band_index][dir_index][cell_index][face_index]
                                        * cell_face_area[cell_index][face_index]
                                        * control_angles[dir_index];
                            }
                            else {
                                heat += ee_curr[cell_index][dir_index][band_index]
                                        * S_dot_normal_cache[band_index][dir_index][cell_index][face_index]
                                        * cell_face_area[cell_index][face_index];
                            }
                        }
                        cell_length += cell_face_area[cell_index][face_index];
                    }
                }
            }
            bc_band_heat_flux[bc_index][band_index] = heat / cell_length * (*bands)[band_index].group_velocity;
            bc_heat_flux[bc_index] += bc_band_heat_flux[bc_index][band_index];
        }
        std::cout << "heat flux of Boundary Condition #" << -1 - bcs->boundaryConditions[bc_index].index << ": " << bc_heat_flux[bc_index] << std::endl;
    }
}

void StaticBTESolver::_preprocess() {
    // get_direction
    if (DM == 2) {
        double delta_phi = 0.5 * PI / num_phi;
        std::vector<double> phi(4 * num_phi, 0);
        phi[0] = 0.5 * delta_phi;
        for (int np = 1; np < phi.size(); np++) {
            phi[np] = phi[np - 1] + delta_phi;
        }
        control_angles.resize(phi.size());
        direction_vectors.reserve(phi.size());
        S.reserve(phi.size());
        for (int np = 0; np < phi.size(); np++) {
            control_angles[np] = WFACTOR / phi.size();
            double x = control_angles[np] * sin(phi[np]);
            double y = control_angles[np] * cos(phi[np]);
            auto sptr = std::make_shared<Point>(x, y, 0);
            S.push_back(sptr);
            x = sin(phi[np]);
            y = cos(phi[np]);
            auto dptr = std::make_shared<Point>(x, y, 0);
            direction_vectors.push_back(dptr);
        }
    }
    else if (DM == 3) {
        double delta_theta = 0.5 * PI / num_theta;
        double delta_phi = 0.5 * PI / num_phi;
        std::vector<double> theta, phi(4 * num_phi, 0);
        if (mesh->dim == 2) {
            theta.resize(num_theta, 0);
        }
        else if (mesh->dim == 3) {
            theta.resize(2 * num_theta, 0);
        }
        else {
            std::cout << "DM is not set properly." << std::endl;
            exit(1);
        }
        theta[0] = 0.5 * delta_theta;
        phi[0] = 0.5 * delta_phi;
        for (int np = 1; np < phi.size(); np++) {
            phi[np] = phi[np - 1] + delta_phi;
        }
        for (int nt = 1; nt < theta.size(); nt++) {
            theta[nt] = theta[nt - 1] + delta_theta;
        }
        control_angles.resize(theta.size() * phi.size());
        direction_vectors.reserve(theta.size() * phi.size());
        S.reserve(theta.size() * phi.size());
        for (int nt = 0; nt < theta.size(); nt++) {
            for (int np = 0; np < phi.size(); np++) {
                int nf = np + nt * 4 * num_phi;
                control_angles[nf] = WFACTOR * 2 * sin(theta[nt]) * sin(0.5 * delta_theta) * delta_phi;
                double x = WFACTOR * sin(phi[np]) * sin(0.5 * delta_phi) *
                           (delta_theta - cos(2 * theta[nt]) * sin(delta_theta));
                double y = WFACTOR * cos(phi[np]) * sin(0.5 * delta_phi) *
                           (delta_theta - cos(2 * theta[nt]) * sin(delta_theta));
                double z = WFACTOR * 0.5 * delta_phi * sin(2 * theta[nt]) * sin(delta_theta);
                auto sptr = std::make_shared<Point>(x, y, z);
                S.push_back(sptr);

                x = sin(theta[nt]) * sin(phi[np]);
                y = sin(theta[nt]) * cos(phi[np]);
                z = cos(theta[nt]);
                auto dptr = std::make_shared<Point>(x, y, z);
                direction_vectors.push_back(dptr);
            }
        }

    }

    // compute cell center and cell volume/area
    cell_centers.reserve(N_cell);
    cell_volume.reserve(N_cell);
    if (mesh->dim == 2) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto p1 = mesh->meshPts[mesh->elements2D[cell_index]->index[0]];
            auto p2 = mesh->meshPts[mesh->elements2D[cell_index]->index[1]];
            auto p3 = mesh->meshPts[mesh->elements2D[cell_index]->index[2]];
            auto ptr = std::make_shared<Point>((p1->x + p2->x + p3->x) / 3,
                                        (p1->y + p2->y + p3->y) / 3);
            cell_centers.push_back(std::move(ptr));
            // FIXME: can be optimized to avoid creating temp object
            // possible fix: implement getArea for shared pointer
            cell_volume.push_back(getArea(*p1, *p2, *p3));
        }
    }
    else if (mesh->dim == 3) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto p1 = mesh->meshPts[mesh->elements3D[cell_index]->index[0]];
            auto p2 = mesh->meshPts[mesh->elements3D[cell_index]->index[1]];
            auto p3 = mesh->meshPts[mesh->elements3D[cell_index]->index[2]];
            auto p4 = mesh->meshPts[mesh->elements3D[cell_index]->index[3]];
            auto ptr = std::make_shared<Point>((p1->x + p2->x + p3->x + p4->x) / 4,
                                        (p1->y + p2->y + p3->y + p4->y) / 4,
                                        (p1->z + p2->z + p3->z + p4->z) / 4);
            cell_centers.push_back(std::move(ptr));
            // FIXME: can be optimized to avoid creating temp object
            // possible fix: implement getVolume for shared pointer
            cell_volume.push_back(getVolume(*p1, *p2, *p3, *p4));
        }
    }
    else {
        std::cout << "DG is not set properly." << std::endl;
        exit(1);
    }

    // compute cell face normal vector and cell face area
    if (mesh->dim == 2) {
        cell_face_normal.resize(N_cell, std::vector<std::shared_ptr<Point>>(N_face));
        cell_face_area.resize(N_cell, std::vector<double>(N_face));
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int edge_index = 0; edge_index < N_face; edge_index++) {
                auto p0 = mesh->meshPts[mesh->elements2D[cell_index]->index[edge_index % 3]];
                auto p1 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 1) % 3]];
                auto p2 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 2) % 3]];
                cell_face_area[cell_index][edge_index] = getLength(p1, p2);
                auto p01 = std::make_shared<Point>(p0->x - p1->x, p0->y - p1->y);
                auto norm = std::make_shared<Point>(p2->y - p1->y, p1->x - p2->x);
                if (dot_prod(norm, p01) > 0) {
                    *norm = Point() - *norm;
                }
                *norm = *norm / norm->length();
                cell_face_normal[cell_index][edge_index] = norm;
            }
        }
    }
    else if (mesh->dim == 3) {
        cell_face_normal.resize(N_cell, std::vector<std::shared_ptr<Point>>(N_face));
        cell_face_area.resize(N_cell, std::vector<double>(N_face));
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int face_index = 0; face_index < N_face; face_index++) {
                auto p0 = mesh->meshPts[mesh->elements3D[cell_index]->index[face_index % 4]];
                auto p1 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 1) % 4]];
                auto p2 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 2) % 4]];
                auto p3 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 3) % 4]];
                cell_face_area[cell_index][face_index] = getArea(p1, p2, p3);
                auto p31 = std::make_shared<Point>();
                auto p32 = std::make_shared<Point>();
                auto p01 = std::make_shared<Point>();
                *p31 = *p3 - *p1;
                *p32 = *p3 - *p2;
                *p01 = *p0 - *p1;
                auto norm = cross_prod(p31, p32);
                if (dot_prod(norm, p01) > 0) {
                    *norm = Point() - *norm;
                }
                *norm = *norm / norm->length();
                cell_face_normal[cell_index][face_index] = norm;
            }
        }
    }
    else {
        std::cout << "DG is not set properly." << std::endl;
        exit(1);
    }

    // compute neighbor cell index
    if (mesh->dim == 2) {
        cell_neighbor_indices.resize(N_cell, std::vector<int>(N_face));
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int edge_index = 0; edge_index < N_face; edge_index++) {
                // CAUTION: order
                int v1 = mesh->elements2D[cell_index]->index[(edge_index + 1) % N_face];
                int v2 = mesh->elements2D[cell_index]->index[(edge_index + 2) % N_face];
                std::unordered_set<int> temp {v1, v2};
                for (int neighbor_index = 0; neighbor_index < N_cell; neighbor_index++) {
                    if (cell_index == neighbor_index) continue;
                    int* neighbor_ptr = mesh->elements2D[neighbor_index]->index;
                    int hit = 0;
                    for (int i = 0; i < N_face; i++) {
                        if (temp.find(neighbor_ptr[i]) != temp.end()) {
                            hit++;
                        }
                    }
                    if (hit >= 2) {
                        cell_neighbor_indices[cell_index][edge_index] = neighbor_index;
                        break;
                    }
                }
                if (cell_neighbor_indices[cell_index][edge_index] == 0) {
                    for (auto& seg_ptr : mesh->elements1D) {
                        int hit = 0;
                        for (int& i : seg_ptr->index) {
                            if (temp.find(i) != temp.end()) {
                                hit++;
                            }
                        }

                        if (hit >= 2) {
                            cell_neighbor_indices[cell_index][edge_index] = - seg_ptr->entity_index - 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    else if (mesh->dim == 3) {
        cell_neighbor_indices.resize(N_cell, std::vector<int>(N_face));
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int face_index = 0; face_index < N_face; face_index++) {
                int v1 = mesh->elements3D[cell_index]->index[(face_index + 1) % N_face];
                int v2 = mesh->elements3D[cell_index]->index[(face_index + 2) % N_face];
                int v3 = mesh->elements3D[cell_index]->index[(face_index + 3) % N_face];
                std::unordered_set<int> temp {v1, v2, v3};
                for (int neighbor_index = 0; neighbor_index < N_cell; neighbor_index++) {
                    if (cell_index == neighbor_index) continue;
                    int* neighbor_ptr = mesh->elements3D[neighbor_index]->index;
                    int hit = 0;
                    for (int i = 0; i < N_face; i++) {
                        if (temp.find(neighbor_ptr[i]) != temp.end()) hit++;
                    }
                    if (hit >= 3) {
                        cell_neighbor_indices[cell_index][face_index] = neighbor_index;
                        break;
                    }
                }
                if (cell_neighbor_indices[cell_index][face_index] == 0) {
                    for (auto& tri_ptr : mesh->elements2D) {
                        int hit = 0;
                        for (int & i : tri_ptr->index) {
                            if (temp.find(i) != temp.end()) {
                                hit++;
                            }
                        }
                        if (hit >= 3) {
                            cell_neighbor_indices[cell_index][face_index] = - tri_ptr->entity_index - 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    else {
        std::cout << "DG is not set properly." << std::endl;
        exit(1);
    }
    _get_const_coefficient();
}

void StaticBTESolver::_iteration(int max_iter) {
    std::cout << "N_cell: " << N_cell << std::endl
              << "N_dir: " << N_dir << std::endl
              << "N_band: " << N_band << std::endl
              << "num_theta: " << num_theta << std::endl
              << "num_phi: " << num_phi << std::endl << std::endl;

    // load inputee or set 0 to start over
    ee_curr = vector3D<double>(N_cell, vector2D<double>(N_dir, std::vector<double>(N_band, 0)));
    for (int iter_index = 0; iter_index < max_iter; iter_index++) {
        ee_prev = ee_curr;
        ee_curr = vector3D<double>(N_cell, vector2D<double>(N_dir, std::vector<double>(N_band, 0)));
        for (int band_index = 0; band_index < N_band; band_index++) {
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                // reconstruct Ke
                vector2D<double> Ke(N_cell, std::vector<double>(N_cell, 0));
                std::vector<double> Re = _get_coefficient(dir_index, band_index);
                for (int i = 0; i < Ke_serialized.size() / 5; i++) {
                    if (Ke_serialized[5 * i] == band_index && Ke_serialized[5 * i + 1] == dir_index) {
                        Ke[Ke_serialized[5 * i + 2]][Ke_serialized[5 * i + 3]] = Ke_serialized[5 * i + 4];
                    }
                }
                std::vector<double> sol = _solve_matrix(Ke, Re);
                for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                    ee_curr[cell_index][dir_index][band_index] = sol[cell_index];
                }
            }
            _get_cell_temperature(band_index);
        }
        _recover_temperature();

        double margin = _get_margin();
        std::cout << "Iteration #" << iter_index << "\t Margin per band per cell: " << margin << std::endl;

        if (margin <= 0.0001) {
            std::cout << std::endl;
            break;
        }
    }
}

void StaticBTESolver::_postprocess() {
    _get_heat_flux();
    std::ofstream outFile;
    outFile.open("Tempcell.dat");
    for (int i = 0; i < N_cell; i++) {
        outFile << (cell_centers[i]->x) / mesh->L_x << " "
                << (cell_centers[i]->y) / mesh->L_y << " ";
        if (mesh->dim == 3) {
            outFile << (cell_centers[i]->z) / mesh->L_z << " ";
        }
        outFile << cell_temperature[i] << std::endl;
    }
    outFile.close();
}

void StaticBTESolver::solve(int max_iter) {
    _preprocess();
    _iteration(max_iter);
    _postprocess();
}