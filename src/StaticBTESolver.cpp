//
// Created by Yucheng Shi on 7/8/20.
//
#include "StaticBTESolver.h"
#include <mpi.h>
#ifdef USE_TIME
#include <chrono>
#endif
#ifdef USE_GPU
#include "bicgstab.hpp"
#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#endif
#else
#include <petscksp.h>
#endif

StaticBTESolver::StaticBTESolver(BTEMesh* mesh, BTEBoundaryCondition* bcs, BTEBand* bands) {
    this->mesh = mesh;
    this->bcs = bcs;
    this->bands = bands;
#ifndef USE_GPU
    PetscInitialize(nullptr, nullptr, (char*)nullptr, nullptr);
#endif
    MPI_Comm_size(MPI_COMM_WORLD, &this->num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->world_rank);
#ifdef USE_GPU
    cudaGetDeviceCount(&this->device_count);
    this->device_id = this->world_rank % this->device_count;
    cudaSetDevice(this->device_id);
    if (this->world_rank == 0) {
        std::cout << this->num_proc << " process(es),   " << this->device_count << " device(s)"<< std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Bind solver rank " << this->world_rank << " to device " << this->device_id << "." << std::endl;
#endif
}


void StaticBTESolver::setParam(int DM, int num_theta, int num_phi, double T_ref) {
    this->num_theta = num_theta;
    this->num_phi = num_phi;
    this->T_ref = T_ref;
    this->DM = DM;
    this->solid_angle = 4 * PI;
    if (mesh->dim == 1) {
        this->T_ref = 0; // need documentation
        N_cell = mesh->elements1D.size();
        if (DM == 3) {
            N_dir = num_theta * num_phi;
        }
        else if (DM == 2) {
            N_dir = num_theta;
            this->solid_angle = 2 * PI;
        }
        N_face = 2;
    }
    if (mesh->dim == 2) {
        N_cell = mesh->elements2D.size();
        if (DM == 3) {
            this->WFACTOR = 2.0;
            N_dir = 4 * num_theta * num_phi;
        }
        else if (DM == 2) {
            this->WFACTOR = 4 * PI;
            N_dir = 4 * num_phi;
        }
        N_face = 3;
    }
    else if (mesh->dim == 3) {
        N_cell = mesh->elements3D.size();
        if (DM == 3) {
            this->WFACTOR = 1.0;
            N_dir = 8 * num_theta * num_phi;
        }
        else if (DM == 2) {
            this->WFACTOR = 2.0;
            N_dir = 4 * num_theta * num_phi;
        }
        N_face = 4;
    }
    N_band = bands->size();
    cell_band_temperature.resize(N_cell, std::vector<double>(N_band));
    cell_band_density.resize(N_cell, std::vector<double>(N_band));
    cell_temperature.resize(N_cell);
}

void StaticBTESolver::_get_cell_temperature(int band_index) {
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        double e0 = 0;
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
            e0 += ee_curr[band_index]->get(dir_index, cell_index) * control_angles[dir_index];
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
            cell_band_density[cell_index][band_index] = (*bands)[band_index].Ctot / (solid_angle) * (cell_temperature[cell_index] - T_ref);
        }
    }
}

void StaticBTESolver::_get_Ke(int band_index, int dir_index, int* csrRowPtr, int* csrColInd, double* csrVal) {
    int csrRowPtr_iter = 0, csrColInd_iter = 0, csrVal_iter = 0;
    csrRowPtr[csrRowPtr_iter++] = 0;
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        std::vector<std::pair<int, double>> compressed_Ke;
        compressed_Ke.reserve(N_face + 1);
        double Ke_cell_index = 0;
        for (int face_index = 0; face_index < N_face; face_index++) {
            double temp = (*bands)[band_index].group_velocity * (*bands)[band_index].relaxation_time *
                          cell_face_area[cell_index][face_index];
            if (mesh->dim > 1) {
                temp *= dot_prod(S[dir_index], cell_face_normal[cell_index][face_index]);
            } else {
                temp *= dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]);
            }
            if (dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]) >= 0) {
                Ke_cell_index += temp;
            } else if (cell_neighbor_indices[cell_index][face_index] >= 0) {
                compressed_Ke.emplace_back(cell_neighbor_indices[cell_index][face_index], temp);
            }
        }
        if (mesh->dim > 1) {
            Ke_cell_index += cell_volume[cell_index] * control_angles[dir_index];
        } else {
            Ke_cell_index += cell_volume[cell_index];
        }

        compressed_Ke.emplace_back(cell_index, Ke_cell_index);

        std::sort(compressed_Ke.begin(), compressed_Ke.end());
        for (auto &p : compressed_Ke) {
            csrColInd[csrColInd_iter++] = p.first;
            csrVal[csrVal_iter++] = p.second;
        }

        csrRowPtr[csrRowPtr_iter] = csrRowPtr[csrRowPtr_iter - 1] + compressed_Ke.size();
        csrRowPtr_iter++;
    }
}

void StaticBTESolver::_get_Re(int band_index, int dir_index, double* Re) {
    for (int cell_index = 0; cell_index < N_cell; cell_index++) {
        for (int face_index = 0; face_index < N_face; face_index++) {
            if (dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]) >= 0) continue;
            double a_f_total = (*bands)[band_index].group_velocity * (*bands)[band_index].relaxation_time * cell_face_area[cell_index][face_index];
            if (mesh->dim > 1) {
                a_f_total *= dot_prod(S[dir_index], cell_face_normal[cell_index][face_index]);
            }
            else {
                a_f_total *= dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]);
            }
            for (auto bc : *bcs) {
                if (cell_neighbor_indices[cell_index][face_index] == bc.index) {
                    if (bc.type == 1) {
                        Re[cell_index] -= (*bands)[band_index].Ctot / (solid_angle)
                                          * (bc.temperature - T_ref)
                                          * a_f_total;

                    } else if (bc.type == 2) {
                        double einsum = 0;
                        double temp = 0;
                        for (int dir_index_inner = 0; dir_index_inner < N_dir; dir_index_inner++) {
                            if (dot_prod(direction_vectors[dir_index_inner], cell_face_normal[cell_index][face_index]) > 0) {
                                if (mesh->dim == 3) {
                                    einsum += ee_prev[band_index]->get(dir_index_inner, cell_index)
                                              * dot_prod(S[dir_index_inner], cell_face_normal[cell_index][face_index])
                                              * control_angles[dir_index_inner];
                                    temp += dot_prod(S[dir_index_inner], cell_face_normal[cell_index][face_index])
                                            * control_angles[dir_index_inner];
                                } else if (mesh->dim == 2) {
                                    einsum += ee_prev[band_index]->get(dir_index_inner, cell_index)
                                              * dot_prod(S[dir_index_inner], cell_face_normal[cell_index][face_index]);
                                    temp += dot_prod(S[dir_index_inner], cell_face_normal[cell_index][face_index]);
                                } else {
                                    einsum += ee_prev[band_index]->get(dir_index_inner, cell_index)
                                              * dot_prod(direction_vectors[dir_index_inner], cell_face_normal[cell_index][face_index])
                                              * control_angles[dir_index_inner];
                                }
                            }
                        }
                        if (mesh->dim != 1) einsum = einsum / temp;
                        else einsum = einsum / PI;
                        Re[cell_index] -= einsum * a_f_total;
                    } else if (bc.type == 31) {
                        int itheta = dir_index / (4 * num_phi);
                        int iphi = dir_index % (4 * num_phi);
                        if (mesh->dim == 2) {
                            if (iphi >= 2 * num_phi) {
                                iphi = 6 * num_phi - iphi - 1;
                            } else {
                                iphi = 2 * num_phi - iphi - 1;
                            }
                        } else if (mesh->dim == 3) {
                            iphi = 4 * num_phi - iphi - 1;
                        }
                        else {
                            std::cout << "Boundary type 31 not supported for DG 1" << std::endl;
                            exit(1);
                        }
                        Re[cell_index] -= ee_prev[band_index]->get(itheta * 4 * num_phi + iphi, cell_index)
                                          * a_f_total;
                    } else if (bc.type == 32) {
                        int itheta = dir_index / (4 * num_phi);
                        int iphi = dir_index % (4 * num_phi);
                        if (mesh->dim == 2) {
                            iphi = 4 * num_phi - iphi - 1;
                            Re[cell_index] -= ee_prev[band_index]->get(itheta * 4 * num_phi + iphi, cell_index)
                                              * a_f_total;
                        } else if (mesh->dim == 3) {
                            int ix = iphi / (2 * num_phi);
                            int isx = iphi % (2 * num_phi);
                            isx = 2 * num_phi - isx - 1;
                            Re[cell_index] -=
                                    ee_prev[band_index]->get(itheta * 4 * num_phi + isx + ix * 2 * num_phi, cell_index)
                                    * a_f_total;
                        }
                        else {
                            std::cout << "Boundary type 31 not supported for DG 1" << std::endl;
                            exit(1);
                        }
                    } else if (bc.type == 33) {
                        if (mesh->dim < 3) {
                            std::cout << "Boundary Condition Type " << bc.type << " not supported. Abort."
                                      << std::endl;
                            exit(1);
                        }
                        int iz = dir_index / (4 * num_phi);
                        int isz = dir_index % (4 * num_phi);
                        iz = 2 * num_phi - iz - 1;
                        Re[cell_index] -= ee_prev[band_index]->get(isz + iz * 4 * num_phi, cell_index)
                                          * a_f_total;
                    } else {
                        std::cout << "Boundary Condition Type " << bc.type << " not supported. Abort." << std::endl;
                        exit(1);
                    }
                }
            }
        }
        if (mesh->dim != 1) {
            Re[cell_index] +=
                    cell_band_density[cell_index][band_index] * cell_volume[cell_index] * control_angles[dir_index];
        }
        else {
            Re[cell_index] +=
                    cell_band_density[cell_index][band_index] * cell_volume[cell_index];
        }
    }
}
#ifndef USE_GPU
double* StaticBTESolver::_solve_matrix(int* RowPtr, int* ColInd, double* Val, double* Re, double* sol) {
    Mat            A;
    Vec            b, x;
    KSP            ksp;
    PC             pc;

    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, N_cell, N_cell, RowPtr, ColInd, Val, &A);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    VecCreateSeq(PETSC_COMM_SELF, N_cell, &b);
    VecDuplicate(b, &x);
    int* indices = new int[N_cell];
    for (int i = 0; i < N_cell; i++) {
        indices[i] = i;
    }
    VecSetValues(b, N_cell, indices, &Re[0], ADD_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
//    PCSetType(pc, PCJACOBI);
//    for mesh dimension 1, it appears ilu preconditioner works better than jacobi pc
//    PCSetType(pc, PCILU);
//    KSPSetType(ksp, KSPBCGS);
    KSPSetTolerances(ksp, 1e-9, PETSC_DEFAULT, PETSC_DEFAULT, 1000);

    KSPSolve(ksp, b, x);

    VecGetValues(x, N_cell, indices, &sol[0]);
    delete [] indices;
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
}
#endif

double StaticBTESolver::_get_margin() {
    double margin = 0;
    for (int band_index = 0; band_index < N_band; band_index++) {
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                margin += std::pow(ee_curr[band_index]->get(dir_index, cell_index) - ee_prev[band_index]->get(dir_index, cell_index), 2);
            }
        }
    }
    return (margin / N_band / N_cell);
}

void StaticBTESolver::_get_heat_flux() {
    bc_band_heat_flux = staticbtesolver::vector2D<double>(bcs->size(), std::vector<double>(N_band, 0));
    bc_heat_flux = std::vector<double>(bcs->size(), 0);
    for (int bc_index = 0; bc_index < bcs->size(); bc_index++) {
        staticbtesolver::BoundaryCondition bc = (*bcs)[bc_index];
        for (int band_index = 0; band_index < N_band; band_index++) {
            double cell_length = 0;
            double heat = 0;
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                for (int face_index = 0; face_index < N_face; face_index++) {
                    if (cell_neighbor_indices[cell_index][face_index] == bc.index) {
                        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                            // TODO: check if control angle is needed
                            if (DM == 3 && mesh->dim == 3) {
                                heat += ee_curr[band_index]->get(dir_index, cell_index)
                                        * dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index])
                                        * cell_face_area[cell_index][face_index]
                                        * control_angles[dir_index];
                            }
                            else if (mesh->dim == 1) {
                                heat += ee_curr[band_index]->get(dir_index, cell_index)
                                        * dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index])
                                        * control_angles[dir_index];
                            }
                            else {
                                heat += ee_curr[band_index]->get(dir_index, cell_index)
                                        * dot_prod(S[dir_index], cell_face_normal[cell_index][face_index])
                                        * cell_face_area[cell_index][face_index];
                            }
                        }
                        cell_length += cell_face_area[cell_index][face_index];
                    }
                }
            }
            if (mesh->dim != 1) {
                bc_band_heat_flux[bc_index][band_index] = heat / cell_length * (*bands)[band_index].group_velocity;
            }
            else {
                bc_band_heat_flux[bc_index][band_index] = heat * (*bands)[band_index].group_velocity;
            }
            bc_heat_flux[bc_index] += bc_band_heat_flux[bc_index][band_index];
        }
        if (this->world_rank == 0) {
            std::cout << "heat flux of Boundary Condition #" << -1 - bcs->boundaryConditions[bc_index].index << ": "
                      << bc_heat_flux[bc_index] << std::endl;
        }
    }
}

void StaticBTESolver::_preprocess() {
    // get_direction
    if (DM == 2) {
        if (mesh->dim >= 2) {
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
                auto sptr = std::make_shared<staticbtesolver::Point>(x, y, 0);
                S.push_back(sptr);
                x = sin(phi[np]);
                y = cos(phi[np]);
                auto dptr = std::make_shared<staticbtesolver::Point>(x, y, 0);
                direction_vectors.push_back(dptr);
            }
        }
        else {
            control_angles.resize(N_dir);
            direction_vectors.reserve(N_dir);
            std::vector<double> cost, sint;
            std::vector<double> W;
            std::vector<double> gauss = staticbtesolver::GaussIntegrationPoints(0, PI, N_dir);;
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                cost.push_back(cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(std::pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                control_angles[dir_index] = W[dir_index] * 2;
                auto dptr = std::make_shared<staticbtesolver::Point>(cost[dir_index], sint[dir_index], 0);
                direction_vectors.push_back(dptr);
            }
            this->WFACTOR = 0;
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                this->WFACTOR += control_angles[dir_index];
            }
        }
    }
    else if (DM == 3) {
        if (mesh->dim >= 2) {
            double delta_theta = 0.5 * PI / num_theta;
            double delta_phi = 0.5 * PI / num_phi;
            std::vector<double> theta, phi(4 * num_phi, 0);
            if (mesh->dim == 2) {
                theta.resize(num_theta, 0);
            } else if (mesh->dim == 3) {
                theta.resize(2 * num_theta, 0);
            } else {
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
                    auto sptr = std::make_shared<staticbtesolver::Point>(x, y, z);
                    S.push_back(sptr);

                    x = sin(theta[nt]) * sin(phi[np]);
                    y = sin(theta[nt]) * cos(phi[np]);
                    z = cos(theta[nt]);
                    auto dptr = std::make_shared<staticbtesolver::Point>(x, y, z);
                    direction_vectors.push_back(dptr);
                }
            }
        }
        else {
            control_angles.resize(N_dir);
            direction_vectors.reserve(N_dir);
            std::vector<double> cost, cosp, sint, sinp;
            std::vector<double> W, Wphi;
            std::vector<double> gauss = staticbtesolver::GaussIntegrationPoints(-1, 1, num_theta);
            std::vector<double> gaussp = staticbtesolver::GaussIntegrationPoints(-PI, PI, num_phi);
            for (int i = 0; i < num_theta; i++) {
                cost.push_back(gauss[2 * i]);
                W.push_back(gauss[2 * i + 1]);
                sint.push_back(std::pow(1 - cost[i] * cost[i], 0.5));
            }
            for (int i = 0; i < num_phi; i++) {
                cosp.push_back(gaussp[2 * i] / PI);
                Wphi.push_back(gaussp[2 * i + 1]);
                sinp.push_back(std::pow(1 - cosp[i] * cosp[i], 0.5));
            }
            for (int i = 0; i< num_theta; i++) {
                for (int j = 0; j < num_phi; j++) {
                    control_angles[i * num_phi + j] = W[i] * Wphi[j];
                    auto dptr = std::make_shared<staticbtesolver::Point>(cost[i], sint[i] * cosp[j], sint[i] * sinp[j]);
                    direction_vectors.push_back(std::move(dptr));
                }
            }
            this->WFACTOR = 0;
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                this->WFACTOR += control_angles[dir_index];
            }
        }
    }

    // compute cell center and cell volume/area
    cell_centers.reserve(N_cell);
    cell_volume.reserve(N_cell);
    std::vector<std::vector<int>> meshPt_2D_map;
    std::vector<std::vector<int>> meshPt_3D_map;
    if (mesh->dim == 1) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto ptr = std::make_shared<staticbtesolver::Point>(mesh->L_x / N_cell * (cell_index + 0.5));
            cell_centers.push_back(std::move(ptr));
            cell_volume.push_back(mesh->L_x / N_cell);
        }
    }
    else if (mesh->dim == 2) {
        meshPt_2D_map.resize(mesh->meshPts.size());
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto p1 = mesh->meshPts[mesh->elements2D[cell_index]->index[0]];
            auto p2 = mesh->meshPts[mesh->elements2D[cell_index]->index[1]];
            auto p3 = mesh->meshPts[mesh->elements2D[cell_index]->index[2]];
            auto ptr = std::make_shared<staticbtesolver::Point>((p1->x + p2->x + p3->x) / 3,
                                        (p1->y + p2->y + p3->y) / 3);
            cell_centers.push_back(std::move(ptr));
            cell_volume.push_back(getArea(*p1, *p2, *p3));
            for (int pt_index : mesh->elements2D[cell_index]->index) {
                meshPt_2D_map[pt_index].push_back(cell_index);
            }
        }
    }
    else if (mesh->dim == 3) {
        meshPt_3D_map.resize(mesh->meshPts.size());
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto p1 = mesh->meshPts[mesh->elements3D[cell_index]->index[0]];
            auto p2 = mesh->meshPts[mesh->elements3D[cell_index]->index[1]];
            auto p3 = mesh->meshPts[mesh->elements3D[cell_index]->index[2]];
            auto p4 = mesh->meshPts[mesh->elements3D[cell_index]->index[3]];
            auto ptr = std::make_shared<staticbtesolver::Point>((p1->x + p2->x + p3->x + p4->x) / 4,
                                        (p1->y + p2->y + p3->y + p4->y) / 4,
                                        (p1->z + p2->z + p3->z + p4->z) / 4);
            cell_centers.push_back(std::move(ptr));
            double temp = getVolume(*p1, *p2, *p3, *p4);
            cell_volume.push_back(temp);
            for (int pt_index : mesh->elements3D[cell_index]->index) {
                meshPt_3D_map[pt_index].push_back(cell_index);
            }
        }
    }
    else {
        std::cout << "DG is not set properly." << std::endl;
        exit(1);
    }

    // compute cell face normal vector and cell face area
    cell_face_normal.resize(N_cell, std::vector<std::shared_ptr<staticbtesolver::Point>>(N_face));
    cell_face_area.resize(N_cell, std::vector<double>(N_face));
    if (mesh->dim == 1) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            cell_face_area[cell_index][0] = 1;
            cell_face_area[cell_index][1] = 1;
            cell_face_normal[cell_index][0] = std::make_shared<staticbtesolver::Point>(-1);
            cell_face_normal[cell_index][1] = std::make_shared<staticbtesolver::Point>(1);
        }
    }
    else if (mesh->dim == 2) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int edge_index = 0; edge_index < N_face; edge_index++) {
                auto p0 = mesh->meshPts[mesh->elements2D[cell_index]->index[edge_index % N_face]];
                auto p1 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 1) % N_face]];
                auto p2 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 2) % N_face]];
                cell_face_area[cell_index][edge_index] = getLength(p1, p2);
                auto p01 = std::make_shared<staticbtesolver::Point>(p0->x - p1->x, p0->y - p1->y);
                auto norm = std::make_shared<staticbtesolver::Point>(p2->y - p1->y, p1->x - p2->x);
                if (dot_prod(norm, p01) > 0) {
                    *norm = staticbtesolver::Point() - *norm;
                }
                *norm = *norm / norm->length();
                cell_face_normal[cell_index][edge_index] = norm;
            }
        }
    }
    else if (mesh->dim == 3) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int face_index = 0; face_index < N_face; face_index++) {
                auto p0 = mesh->meshPts[mesh->elements3D[cell_index]->index[face_index % 4]];
                auto p1 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 1) % 4]];
                auto p2 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 2) % 4]];
                auto p3 = mesh->meshPts[mesh->elements3D[cell_index]->index[(face_index + 3) % 4]];
                cell_face_area[cell_index][face_index] = getArea(p1, p2, p3);
                auto p31 = std::make_shared<staticbtesolver::Point>();
                auto p32 = std::make_shared<staticbtesolver::Point>();
                auto p01 = std::make_shared<staticbtesolver::Point>();
                *p31 = *p3 - *p1;
                *p32 = *p3 - *p2;
                *p01 = *p0 - *p1;
                auto norm = cross_prod(p31, p32);
                if (dot_prod(norm, p01) > 0) {
                    *norm = staticbtesolver::Point() - *norm;
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
    cell_neighbor_indices.resize(N_cell, std::vector<int>(N_face));
    if (mesh->dim == 1) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            cell_neighbor_indices[cell_index][0] = cell_index - 1;
            cell_neighbor_indices[cell_index][1] = cell_index + 1;
        }
        // left boundary hard coded to boundary condition 0
        // right boundary hard coded to boundary condition 1
        cell_neighbor_indices[0][0] = - 0 - 1;
        cell_neighbor_indices[N_cell - 1][1] = - 1 - 1;
    }
    else if (mesh->dim == 2) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int edge_index = 0; edge_index < N_face; edge_index++) {
                int v1 = mesh->elements2D[cell_index]->index[(edge_index + 1) % N_face];
                int v2 = mesh->elements2D[cell_index]->index[(edge_index + 2) % N_face];
                for (int neighbor_index : meshPt_2D_map[v1]) {
                    if (cell_index == neighbor_index) continue;
                    int* neighbor_ptr = mesh->elements2D[neighbor_index]->index;
                    int hit = 0;
                    for (int i = 0; i < N_face; i++) {
                        if (neighbor_ptr[i] == v1 || neighbor_ptr[i] == v2) {
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
                            if (i == v1 || i == v2) {
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
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int face_index = 0; face_index < N_face; face_index++) {
                int v1 = mesh->elements3D[cell_index]->index[(face_index + 1) % N_face];
                int v2 = mesh->elements3D[cell_index]->index[(face_index + 2) % N_face];
                int v3 = mesh->elements3D[cell_index]->index[(face_index + 3) % N_face];
                for (int neighbor_index : meshPt_3D_map[v1]) {
                    if (cell_index == neighbor_index) continue;
                    int* neighbor_ptr = mesh->elements3D[neighbor_index]->index;
                    int hit = 0;
                    for (int i = 0; i < N_face; i++) {
                        if (neighbor_ptr[i] == v1 || neighbor_ptr[i] == v2 || neighbor_ptr[i] == v3) hit++;
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
                            if (i == v1 || i == v2 || i == v3) {
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
}

void StaticBTESolver::_iteration(int max_iter) {
    if (this->world_rank == 0) {
        std::cout << "N_cell: " << N_cell << std::endl
                  << "N_dir: " << N_dir << std::endl
                  << "N_band: " << N_band << std::endl
                  << "num_theta: " << num_theta << std::endl
                  << "num_phi: " << num_phi << std::endl << std::endl;
    }
#ifdef USE_GPU
    auto bicgstab_solver = new BICGSTAB(N_cell, 1000, 1e-9);
#endif

    ee_curr.resize(N_band, nullptr);
    for (int band_index = 0; band_index < N_band; band_index++) {
        ee_curr[band_index] = new staticbtesolver::ContinuousArray(N_dir, N_cell);
    }
    ee_prev.resize(N_band, nullptr);
    for (int iter_index = 0; iter_index < max_iter; iter_index++) {
#ifdef USE_TIME
        auto time_start = std::chrono::high_resolution_clock::now();
#endif
        for (int band_index = 0; band_index < N_band; band_index++) {
            delete ee_prev[band_index];
        }
        ee_prev = ee_curr;
        for (int band_index = 0; band_index < N_band; band_index++) {
            ee_curr[band_index] = new staticbtesolver::ContinuousArray(N_dir, N_cell);
        }
#ifdef USE_TIME
        auto get_Re_time = std::chrono::microseconds(0);
        auto get_Ke_time = std::chrono::microseconds(0);
        auto solver_time = std::chrono::microseconds(0);
        auto start = std::chrono::high_resolution_clock::now();
#endif
        for (int band_index = 0; band_index < N_band; band_index++) {
            for (int dir_index = this->world_rank; dir_index < N_dir; dir_index += this->num_proc) {

                auto csrRowPtr = new int[N_cell + 1];
                auto csrColInd = new int[N_face * N_cell + 1];
                auto csrVal = new double[N_face * N_cell + 1];
#ifdef USE_TIME
                auto _get_Ke_start = std::chrono::high_resolution_clock::now();
#endif
                _get_Ke(band_index, dir_index, csrRowPtr, csrColInd, csrVal);
#ifdef USE_TIME
                auto _get_Ke_end = std::chrono::high_resolution_clock::now();
                get_Ke_time += std::chrono::duration_cast<std::chrono::microseconds>(_get_Ke_end - _get_Ke_start);
#endif
                auto* Re = new double[N_cell]();
                _get_Re(band_index, dir_index, Re);
#ifdef USE_TIME
                auto _get_Re_end = std::chrono::high_resolution_clock::now();
                get_Re_time += std::chrono::duration_cast<std::chrono::microseconds>(_get_Re_end - _get_Ke_end);
                auto solver_start = std::chrono::high_resolution_clock::now();
#endif
                auto* sol = new double[N_cell];
#ifdef USE_GPU
                int nnz = csrRowPtr[N_cell];
                bicgstab_solver->init(csrRowPtr, csrColInd, csrVal, nnz, Re);
                bicgstab_solver->solve(sol);
#else
                _solve_matrix(csrRowPtr, csrColInd, csrVal, Re, sol);
#endif
#ifdef USE_TIME
                auto solver_end = std::chrono::high_resolution_clock::now();
//                std::cout << "[" << band_index << "] Process " << this->world_rank << " is assigned to dir_index = " << dir_index << std::endl;
//                std::cout << "     takes " << std::chrono::duration_cast<std::chrono::milliseconds>(solver_end - solver_start).count() << "ms to solve the system" << std::endl;
                solver_time += std::chrono::duration_cast<std::chrono::microseconds>(solver_end - solver_start);
                auto gather_start = std::chrono::high_resolution_clock::now();
#endif
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Allgather(sol,
                              N_cell,
                              MPI_DOUBLE,
                              (ee_curr[band_index]->data + N_cell * (dir_index - this->world_rank)),
                              N_cell,
                              MPI_DOUBLE,
                              MPI_COMM_WORLD
                              );
#ifdef USE_TIME
                auto gather_end = std::chrono::high_resolution_clock::now();
//                std::cout << "[" << band_index << "] Process " << this->world_rank << " takes " << std::chrono::duration_cast<std::chrono::milliseconds>(gather_end - gather_start).count() << "ms to gather" << std::endl;
#endif
                delete [] sol;
                delete [] Re;
                delete [] csrRowPtr;
                delete [] csrColInd;
                delete [] csrVal;
            }
            _get_cell_temperature(band_index);
        }
#ifdef USE_TIME
        auto time_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);
#endif
        _recover_temperature();

        double margin = _get_margin();
        if (this->world_rank == 0) {
            std::cout << "----------------------------------------------------------------------------------" << std::endl;
            std::cout << "Iteration #" << iter_index << "\t Margin per band per cell: " << margin << std::endl;
#ifdef USE_TIME
            std::cout << "Time taken by inner loop: " << 1.0 * duration.count() / 1000 << " milliseconds" << std::endl;
            std::cout << "Time taken by get_Ke: " << 1.0 * get_Ke_time.count() / 1000 << " milliseconds" << std::endl;
            std::cout << "Time taken by get_Re: " << 1.0 * get_Re_time.count() / 1000 << " milliseconds" << std::endl;
            std::cout << "Time taken by iterative solver: " << 1.0 * solver_time.count() / 1000 << " milliseconds" << std::endl;
#endif
        }
        if (margin <= 0.00001) {
            if (this->world_rank == 0) {
                std::cout << std::endl;
            }
            break;
        }
    }
#ifdef USE_GPU
    delete bicgstab_solver;
#endif
}

void StaticBTESolver::_postprocess() {
    MPI_Barrier(MPI_COMM_WORLD);
    _get_heat_flux();
    MPI_Barrier(MPI_COMM_WORLD);
    for (int band_index = 0; band_index < N_band; band_index++) {
        delete ee_curr[band_index];
        delete ee_prev[band_index];
    }
//if (this->world_rank == 0) {
//    std::cout << std::endl;
//    std::ofstream outFile;
//    outFile.open("Tempcell.dat");
//   for (int i = 0; i < N_cell; i++) {
//      outFile << (cell_centers[i]->x) / mesh->L_x << " "
//              << (cell_centers[i]->y) / mesh->L_y << " ";
//      if (mesh->dim == 3) {
//          outFile << (cell_centers[i]->z) / mesh->L_z << " ";
//      }
//      outFile << cell_temperature[i] << std::endl;
//  }
//   outFile.close();
//}
#ifndef USE_GPU
    PetscFinalize();
#endif
}

#ifdef USE_GPU
size_t StaticBTESolver::print_host_mem() {
    size_t host_mem_before = staticbtesolver::get_host_memory();
    size_t *mem_sum = nullptr;
    if (this->world_rank == 0) {
        mem_sum = new size_t[this->num_proc];
    }
    MPI_Gather(&host_mem_before, 1, my_MPI_SIZE_T, mem_sum, 1, my_MPI_SIZE_T, 0,
               MPI_COMM_WORLD);
    size_t mem = 0;
    if (this->world_rank == 0) {
        for (int i = 0; i < this->num_proc; i++) {
            mem += mem_sum[i];
        }
        delete [] mem_sum;
    }
    return mem / 1024;
}

size_t StaticBTESolver::print_device_mem() {
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    return (total - free) / 1024;
}
#endif

void StaticBTESolver::solve(int max_iter) {
    _preprocess();
#ifdef USE_TIME
    auto start = std::chrono::high_resolution_clock::now();
#endif
    _iteration(max_iter);
#ifdef USE_TIME
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Time taken by iteration: "
         << duration.count() * 0.001 << " milliseconds" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    _postprocess();
}

