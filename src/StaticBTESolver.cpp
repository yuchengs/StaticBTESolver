//
// Created by Yucheng Shi on 7/8/20.
//
#define USE_TIME
#include "StaticBTESolver.h"
#ifdef USE_TIME
#include <chrono>
#endif
#ifdef USE_GPU
#define VIENNACL_WITH_CUDA 1
#include <mpi.h>
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
#else
#include "petscksp.h"
#endif


StaticBTESolver::StaticBTESolver(BTEMesh* mesh, BTEBoundaryCondition* bcs, BTEBand* bands) {
    this->mesh = mesh;
    this->bcs = bcs;
    this->bands = bands;

#ifdef USE_GPU
    MPI_Comm_size(MPI_COMM_WORLD, &this->num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->world_rank);

//    if (mesh->dim == 1) {
//        if (this->world_rank == 0) {
//            std::cout << "GPU version does not support 1D geometry." << std::endl;
//        }
//        exit(0);
//    }

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
            e0 += ee_curr[band_index].get(dir_index, cell_index) * control_angles[dir_index];
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

void StaticBTESolver::_get_const_coefficient() {
    dv_dot_normal_cache.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));
    if (mesh->dim > 1) {
        S_dot_normal_cache.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));
    }
    a_f_total.resize(N_band, vector3D<double>(N_dir, vector2D<double>(N_cell, std::vector<double>())));

#ifdef USE_GPU
    csrRowPtr.resize(N_band, std::vector<unsigned int*>(N_dir, nullptr));
    csrColInd.resize(N_band, std::vector<unsigned int*>(N_dir, nullptr));
    csrVal.resize(N_band, std::vector<double*>(N_dir, nullptr));
#else
    Ke_serialized.resize(N_band, vector2D<double>(N_dir, std::vector<double>()));
#endif

    for (int band_index = 0; band_index < N_band; band_index++) {
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
#ifdef USE_GPU
            std::vector<unsigned int> csrRowPtr_stl;
            csrRowPtr_stl.push_back(0);
            std::vector<unsigned int> csrColInd_stl;
            std::vector<double> csrVal_stl;
#endif
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                std::vector<double> Ke(N_cell, 0);
                for (int face_index = 0; face_index < N_face; face_index++) {
                    double area = cell_face_area[cell_index][face_index];
                    dv_dot_normal_cache[band_index][dir_index][cell_index].push_back(
                            dot_prod(direction_vectors[dir_index], cell_face_normal[cell_index][face_index]));
                    if (mesh->dim > 1) {
                        S_dot_normal_cache[band_index][dir_index][cell_index].push_back(dot_prod(S[dir_index], cell_face_normal[cell_index][face_index]));
                    }
                    double temp = (*bands)[band_index].group_velocity * (*bands)[band_index].relaxation_time * area;
                    if (mesh->dim > 1) {
                        temp *= S_dot_normal_cache[band_index][dir_index][cell_index][face_index];
                    }
                    else {
                        temp *= dv_dot_normal_cache[band_index][dir_index][cell_index][face_index];
                    }
                    if (dv_dot_normal_cache[band_index][dir_index][cell_index][face_index] >= 0) {
                        Ke[cell_index] += temp;
                    } else if (cell_neighbor_indices[cell_index][face_index] >= 0) {
                        Ke[cell_neighbor_indices[cell_index][face_index]] += temp;
                    }
                    a_f_total[band_index][dir_index][cell_index].push_back(temp);
                }
                if (mesh->dim > 1) {
                    Ke[cell_index] += cell_volume[cell_index] * control_angles[dir_index];
                }
                else {
                    Ke[cell_index] += cell_volume[cell_index];
                }
#ifdef USE_GPU
                std::vector<int> temp_indices;
                for (int face_index = 0; face_index < N_face; face_index++) {
                    if (cell_neighbor_indices[cell_index][face_index] >= 0 && Ke[cell_neighbor_indices[cell_index][face_index]] != 0) {
                        temp_indices.push_back(cell_neighbor_indices[cell_index][face_index]);
                    }
                }
                if (Ke[cell_index] != 0) temp_indices.push_back(cell_index);
                std::sort(temp_indices.begin(), temp_indices.end());

                for (int index : temp_indices) {
                    csrColInd_stl.push_back(index);
                    csrVal_stl.push_back(Ke[index]);
                }
                csrRowPtr_stl.push_back(csrRowPtr_stl.back() + temp_indices.size());
#else
                if (Ke[cell_index] != 0) {
                    Ke_serialized[band_index][dir_index].push_back(cell_index);
                    Ke_serialized[band_index][dir_index].push_back(cell_index);
                    Ke_serialized[band_index][dir_index].push_back(Ke[cell_index]);
                }
                for (int face_index = 0; face_index < N_face; face_index++) {
                    if (cell_neighbor_indices[cell_index][face_index] >= 0) {
                        Ke_serialized[band_index][dir_index].push_back(cell_index);
                        Ke_serialized[band_index][dir_index].push_back(cell_neighbor_indices[cell_index][face_index]);
                        Ke_serialized[band_index][dir_index].push_back(Ke[cell_neighbor_indices[cell_index][face_index]]);
                    }
                }
#endif
            }
#ifdef USE_GPU
            auto* csrRowPtr_c = new unsigned int[csrRowPtr_stl.size()];
            std::copy(csrRowPtr_stl.begin(), csrRowPtr_stl.end(), csrRowPtr_c);
            auto* csrColInd_c = new unsigned int[csrColInd_stl.size()];
            std::copy(csrColInd_stl.begin(), csrColInd_stl.end(), csrColInd_c);
            auto* csrVal_c = new double[csrVal_stl.size()];
            std::copy(csrVal_stl.begin(), csrVal_stl.end(), csrVal_c);
            csrRowPtr[band_index][dir_index] = csrRowPtr_c;
            csrColInd[band_index][dir_index] = csrColInd_c;
            csrVal[band_index][dir_index] = csrVal_c;
//            std::cout << "nnz = " << " " << csrRowPtr_c[N_cell] << std::endl;
//            std::cout << "ColInd: " << std::endl;
//            for (int i = 0; i < csrRowPtr_c[N_cell]; i++) {
//                std::cout << csrColInd_c[i] <<  " ";
//            }
//            std::cout << std::endl << "Val:" << std::endl;
//            for (int i = 0; i < csrRowPtr_c[N_cell]; i++) {
//                std::cout << csrVal_c[i] <<  " ";
//            }
//            std::cout << std::endl << "RowPtr:" << std::endl;
//            for (int i = 0; i < N_cell; i++) {
//                std::cout << csrRowPtr_c[i] <<  " ";
//            }
//            std::cout << std::endl;
//            exit(0);
#endif
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
                        Re[cell_index] -= (*bands)[band_index].Ctot / (solid_angle)
                                          * (bc.temperature - T_ref)
                                          * a_f_total[band_index][dir_index][cell_index][face_index];

                    } else if (bc.type == 2) {
                        double einsum = 0;
                        double temp = 0;
                        for (int dir_index_inner = 0; dir_index_inner < N_dir; dir_index_inner++) {
                            if (dv_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index] > 0) {
                                if (mesh->dim == 3) {
                                    einsum += ee_prev[band_index].get(dir_index, cell_index)
                                              * S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index]
                                              * control_angles[dir_index_inner];
                                    temp += S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index]
                                            * control_angles[dir_index_inner];
                                } else if (mesh->dim == 2) {
                                    einsum += ee_prev[band_index].get(dir_index, cell_index)
                                              * S_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index];
                                } else {
                                    einsum += ee_prev[band_index].get(dir_index, cell_index)
                                              * dv_dot_normal_cache[band_index][dir_index_inner][cell_index][face_index]
                                              * control_angles[dir_index_inner];
                                }
                            }
                        }
                        if (mesh->dim == 3) einsum = einsum / temp;
                        else einsum = einsum / PI;
                        Re[cell_index] -= einsum * a_f_total[band_index][dir_index][cell_index][face_index];
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
                        Re[cell_index] -= ee_prev[band_index].get(itheta * 4 * num_phi + iphi, cell_index)
                                          * a_f_total[band_index][dir_index][cell_index][face_index];
                    } else if (bc.type == 32) {
                        int itheta = dir_index / (4 * num_phi);
                        int iphi = dir_index % (4 * num_phi);
                        if (mesh->dim == 2) {
                            iphi = 4 * num_phi - iphi - 1;
                            Re[cell_index] -= ee_prev[band_index].get(itheta * 4 * num_phi + iphi, cell_index)
                                              * a_f_total[band_index][dir_index][cell_index][face_index];
                        } else if (mesh->dim == 3) {
                            int ix = iphi / (2 * num_phi);
                            int isx = iphi % (2 * num_phi);
                            isx = 2 * num_phi - isx - 1;
                            Re[cell_index] -=
                                    ee_prev[band_index].get(itheta * 4 * num_phi + isx + ix * 2 * num_phi, cell_index)
                                    * a_f_total[band_index][dir_index][cell_index][face_index];
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
                        Re[cell_index] -= ee_prev[band_index].get(isz + iz * 4 * num_phi, cell_index)
                                          * a_f_total[band_index][dir_index][cell_index][face_index];
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
    return Re;
}
#ifndef USE_GPU
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
#endif

double StaticBTESolver::_get_margin() {
    double margin = 0;
    for (int band_index = 0; band_index < N_band; band_index++) {
        for (int dir_index = 0; dir_index < N_dir; dir_index++) {
            for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                margin += std::pow(ee_curr[band_index].get(dir_index, cell_index) - ee_prev[band_index].get(dir_index, cell_index), 2);
            }
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
                                heat += ee_curr[band_index].get(dir_index, cell_index)
                                        * dv_dot_normal_cache[band_index][dir_index][cell_index][face_index]
                                        * cell_face_area[cell_index][face_index]
                                        * control_angles[dir_index];
                            }
                            else if (mesh->dim == 1) {
                                heat += ee_curr[band_index].get(dir_index, cell_index)
                                        * dv_dot_normal_cache[band_index][dir_index][cell_index][face_index]
                                        * control_angles[dir_index];
                            }
                            else {
                                heat += ee_curr[band_index].get(dir_index, cell_index)
                                        * S_dot_normal_cache[band_index][dir_index][cell_index][face_index]
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
#ifdef USE_GPU
        if (this->world_rank == 0) {
#endif
            std::cout << "heat flux of Boundary Condition #" << -1 - bcs->boundaryConditions[bc_index].index << ": "
                      << bc_heat_flux[bc_index] << std::endl;
#ifdef USE_GPU
        }
#endif
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
                auto sptr = std::make_shared<Point>(x, y, 0);
                S.push_back(sptr);
                x = sin(phi[np]);
                y = cos(phi[np]);
                auto dptr = std::make_shared<Point>(x, y, 0);
                direction_vectors.push_back(dptr);
            }
        }
        else {
            control_angles.resize(N_dir);
            direction_vectors.reserve(N_dir);
            std::vector<double> cost, sint;
            std::vector<double> W;
            std::vector<double> gauss = GaussIntegrationPoints(0, PI, N_dir);;
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                cost.push_back(cos(gauss[2 * dir_index]));
                W.push_back(gauss[2 * dir_index + 1]);
                sint.push_back(std::pow(1 - cost[dir_index] * cost[dir_index], 0.5));
            }
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
                control_angles[dir_index] = W[dir_index] * 2;
                auto dptr = std::make_shared<Point>(cost[dir_index], sint[dir_index], 0);
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
        else {
            control_angles.resize(N_dir);
            direction_vectors.reserve(N_dir);
            std::vector<double> cost, cosp, sint, sinp;
            std::vector<double> W, Wphi;
            std::vector<double> gauss = GaussIntegrationPoints(-1, 1, num_theta);
            std::vector<double> gaussp = GaussIntegrationPoints(-PI, PI, num_phi);
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
                    auto dptr = std::make_shared<Point>(cost[i], sint[i] * cosp[j], sint[i] * sinp[j]);
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
    if (mesh->dim == 1) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            auto ptr = std::make_shared<Point>(mesh->L_x / N_cell * (cell_index + 0.5));
            cell_centers.push_back(std::move(ptr));
            cell_volume.push_back(mesh->L_x / N_cell);
        }
    }
    else if (mesh->dim == 2) {
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
            double temp = getVolume(*p1, *p2, *p3, *p4);
            cell_volume.push_back(temp);
        }
    }
    else {
        std::cout << "DG is not set properly." << std::endl;
        exit(1);
    }

    // compute cell face normal vector and cell face area
    cell_face_normal.resize(N_cell, std::vector<std::shared_ptr<Point>>(N_face));
    cell_face_area.resize(N_cell, std::vector<double>(N_face));
    if (mesh->dim == 1) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            cell_face_area[cell_index][0] = 1;
            cell_face_area[cell_index][1] = 1;
            cell_face_normal[cell_index][0] = std::make_shared<Point>(-1);
            cell_face_normal[cell_index][1] = std::make_shared<Point>(1);
        }
    }
    else if (mesh->dim == 2) {
        for (int cell_index = 0; cell_index < N_cell; cell_index++) {
            for (int edge_index = 0; edge_index < N_face; edge_index++) {
                auto p0 = mesh->meshPts[mesh->elements2D[cell_index]->index[edge_index % N_face]];
                auto p1 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 1) % N_face]];
                auto p2 = mesh->meshPts[mesh->elements2D[cell_index]->index[(edge_index + 2) % N_face]];
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
                for (int neighbor_index = 0; neighbor_index < N_cell; neighbor_index++) {
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
                for (int neighbor_index = 0; neighbor_index < N_cell; neighbor_index++) {
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
    _get_const_coefficient();
}

void StaticBTESolver::_iteration(int max_iter) {
#ifdef USE_GPU
    if (this->world_rank == 0) {
#endif
        std::cout << "N_cell: " << N_cell << std::endl
                  << "N_dir: " << N_dir << std::endl
                  << "N_band: " << N_band << std::endl
                  << "num_theta: " << num_theta << std::endl
                  << "num_phi: " << num_phi << std::endl << std::endl;
#ifdef USE_GPU
    }


    viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
    chow_patel_ilu_config.sweeps(3);       // three nonlinear sweeps
    chow_patel_ilu_config.jacobi_iters(2); // two Jacobi iterations per triangular 'solve' Rx=r
#endif

    ee_curr.resize(N_band, ContinuousArray(N_dir, N_band));
    for (int iter_index = 0; iter_index < max_iter; iter_index++) {
#ifdef USE_TIME
        auto time_start = std::chrono::high_resolution_clock::now();
#endif
        ee_prev = ee_curr;
        for (int band_index = 0; band_index < N_band; band_index++) {
            ee_curr[band_index].init(N_dir, N_cell);
        }
        // TIMEDEBUG:
        auto start = std::chrono::high_resolution_clock::now();
        for (int band_index = 0; band_index < N_band; band_index++) {
#ifdef USE_GPU
            for (int dir_index = this->world_rank; dir_index < N_dir; dir_index += this->num_proc) {
#else
            for (int dir_index = 0; dir_index < N_dir; dir_index++) {
#endif
                std::vector<double> Re = _get_coefficient(dir_index, band_index);
#ifdef USE_GPU
                unsigned int nnz = csrRowPtr[band_index][dir_index][N_cell];
                viennacl::compressed_matrix<double> vKe;
                vKe.set(csrRowPtr[band_index][dir_index], csrColInd[band_index][dir_index], csrVal[band_index][dir_index], N_cell, N_cell, nnz);
                viennacl::vector<double> vRe(Re.size());
                viennacl::copy(Re, vRe);

                viennacl::vector<double> vsol;

                if (mesh->dim == 1) {
                    vsol = viennacl::linalg::solve(vKe, vRe, viennacl::linalg::gmres_tag(1e-50, 10000, 30));
                }
                else {
                    viennacl::linalg::chow_patel_ilu_precond<viennacl::compressed_matrix<double>> chow_patel_ilu(vKe, chow_patel_ilu_config);
                    vsol = viennacl::linalg::solve(vKe, vRe,viennacl::linalg::bicgstab_tag(1e-9, 1000,200), chow_patel_ilu);
                }

                auto* sol = new double[N_cell];
                viennacl::copy(vsol.begin(), vsol.end(), sol);
//                for (int i = 0; i < N_cell; i++){
//                    std::cout << sol[i] << " ";
//                }
//                exit(0);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Allgather(sol,
                              N_cell,
                              MPI_DOUBLE,
                              (ee_curr[band_index].data + N_cell * (dir_index - this->world_rank)),
                              N_cell,
                              MPI_DOUBLE,
                              MPI_COMM_WORLD
                              );
                //ee_curr[band_index].print();
                delete [] sol;
#else
                // reconstruct Ke
                vector2D<double> Ke(N_cell, std::vector<double>(N_cell, 0));
                // TODO #1 (possible optimization) Ke_serialized
                for (int i = 0; i < Ke_serialized[band_index][dir_index].size() / 3; i++) {
                    Ke[Ke_serialized[band_index][dir_index][3 * i]][Ke_serialized[band_index][dir_index][3 * i + 1]] = Ke_serialized[band_index][dir_index][3 * i + 2];
                }
//                for (int i = 0; i < N_cell; i++) {
//                    for (int j = 0; j < N_cell; j++) {
//                        std::cout << Ke[i][j] << " ";
//                    }
//                    std::cout << std::endl;
//                }
//                exit(0);
//
//                for (int i = 0; i < N_cell; i++) {
//                    std::cout << sol[i] << " ";
//                }
//                exit(0);
                std::vector<double> sol = _solve_matrix(Ke, Re);
                for (int cell_index = 0; cell_index < N_cell; cell_index++) {
                    *ee_curr[band_index].get_ptr(dir_index, cell_index) = sol[cell_index];
                }
#endif
            }
            _get_cell_temperature(band_index);
        }
        _recover_temperature();

        double margin = _get_margin();
#ifdef USE_TIME
        auto time_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);
#endif
#ifdef USE_GPU
        if (this->world_rank == 0) {
#endif
            std::cout << "Iteration #" << iter_index << "\t Margin per band per cell: " << margin << std::endl;
#ifdef USE_TIME
            std::cout << "Time taken by inner loop: " << 1.0 * duration.count() / 1000 << " milliseconds" << std::endl;
#endif
#ifdef USE_GPU
        }
#endif
        if (margin <= 0.0001) {
#ifdef USE_GPU
            if (this->world_rank == 0) {
#endif
            std::cout << std::endl;
#ifdef USE_GPU
            }
#endif
            break;
        }
    }
}

void StaticBTESolver::_postprocess() {
    _get_heat_flux();
#ifdef USE_GPU
    if (this->world_rank == 0) {
#endif
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
#ifdef USE_GPU
    }
#endif
/*    for (int i = 0; i < N_band; i++) {
        for (int j = 0; j < N_dir; j++) {
            delete [] csrRowPtr[i][j];
            delete [] csrColInd[i][j];
            delete [] csrVal[i][j];
        }
    } */
}

void StaticBTESolver::solve(int max_iter) {
    _preprocess();
#ifdef USE_TIME
    auto start = std::chrono::high_resolution_clock::now();
#endif
    _iteration(max_iter);
#ifdef USE_TIME
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << std::endl << "Time taken by iteration: "
         << duration.count() * 0.001 << " milliseconds" << std::endl;
#endif
    _postprocess();
}

