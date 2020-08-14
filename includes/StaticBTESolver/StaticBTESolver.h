//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_STATICBTESOLVER_H
#define BTESOLVER_STATICBTESOLVER_H

#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include "StaticBTESolver/BTEMesh.h"
#include "StaticBTESolver/BTEBand.h"
#include "StaticBTESolver/BTEBoundaryCondition.h"

class StaticBTESolver {
    int num_proc;
    int world_rank;
    int device_count;
    int device_id;
#ifdef USE_GPU
    int* dev_num_proc_ptr;
    int* dev_world_rank_ptr;
    int* dev_device_count_ptr;
    int* device_id_ptr;
#endif
    // cache
    BTEMesh* mesh;
    BTEBoundaryCondition* bcs;
    BTEBand* bands;
    // CADOM param
    int num_theta;
    int num_phi;
    double WFACTOR;
#ifdef USE_GPU
    int* dev_num_theta_ptr;
    int* dev_num_phi_ptr;
    double* dev_WFACTOR_ptr;
#endif
    // hyper
    int DM;
    double T_ref;
    // intermediate variables
    int N_cell, N_dir, N_band, N_face;
#ifdef USE_GPU
    int* dev_DM_ptr;
    int* dev_T_ref_ptr;
    int* dev_N_cell_ptr;
    int* dev_N_dir_ptr;
    int* dev_N_band_ptr;
    int* dev_N_face_ptr;
#endif
    double solid_angle;
    std::vector<double> control_angles;
    std::vector<std::shared_ptr<Point>> S;
    std::vector<std::shared_ptr<Point>> direction_vectors;
    std::vector<std::shared_ptr<Point>> cell_centers;
    std::vector<double> cell_volume;
    std::vector<std::vector<std::shared_ptr<Point>>> cell_face_normal;
    std::vector<std::vector<double>> cell_face_area;
    std::vector<std::vector<int>> cell_neighbor_indices;
    vector2D<double> cell_band_temperature;
    vector2D<double> cell_band_density;
    std::vector<double> cell_temperature;
    vector2D<unsigned int*> csrRowPtr;
    vector2D<unsigned int*> csrColInd;
    vector2D<double*> csrVal;

    std::vector<ContinuousArray*> ee_curr, ee_prev;
    vector2D<double> bc_band_heat_flux;
    std::vector<double> bc_heat_flux;

    // some utility functions
    void _get_cell_temperature(int band_index);
    void _recover_temperature();
    void _get_const_coefficient();
    std::vector<double> _get_coefficient(int dir_index, int band_index);

    double _get_margin();
    void _get_heat_flux();

    void _preprocess();
    void _iteration(int max_iter);
    void _postprocess();

#ifndef USE_GPU
    std::vector<double> _solve_matrix(int* csrRowPtr, int* csrColInd, double* csrVal, std::vector<double>& Re);
#else
    size_t print_host_mem();
    size_t print_device_mem();
#endif

public:
    StaticBTESolver(BTEMesh* mesh, BTEBoundaryCondition* bcs, BTEBand* bands);
    void setParam(int DM, int num_theta, int num_phi, double T_ref);
    void solve(int max_iter);
};

#endif //BTESOLVER_STATICBTESOLVER_H
