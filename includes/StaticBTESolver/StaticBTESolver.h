//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_STATICBTESOLVER_H
#define BTESOLVER_STATICBTESOLVER_H

#include <vector>
#include <memory>
#include "StaticBTESolver/BTEMesh.h"
#include "StaticBTESolver/BTEBand.h"
#include "StaticBTESolver/BTEBoundaryCondition.h"
#include "utility.h"

class StaticBTESolver {
    int num_proc;
    int world_rank;
    int device_count;
    int device_id;

    // cache
    BTEMesh* mesh;
    BTEBoundaryCondition* bcs;
    BTEBand* bands;
    // CADOM param
    int num_theta;
    int num_phi;
    double WFACTOR;

    // hyper
    int DM;
    double T_ref;
    // intermediate variables
    int N_cell, N_dir, N_band, N_face;

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

    std::vector<ContinuousArray*> ee_curr, ee_prev;
    vector2D<double> bc_band_heat_flux;
    std::vector<double> bc_heat_flux;

    // some utility functions
    void _get_Re(int band_index, int dir_index, double* Re); // N_cell * N_face
    void _get_Ke(int band_index, int dir_index, int* csrRowPtr, int* csrColInd, double* csrVal); // N_cell * N_face
    void _get_cell_temperature(int band_index);
    void _recover_temperature();
    double _get_margin();
    void _get_heat_flux();

    void _preprocess();
    void _iteration(int max_iter);
    void _postprocess();

#ifndef USE_GPU
    double* _solve_matrix(int* csrRowPtr, int* csrColInd, double* csrVal, double* Re, double* sol);
#else
    size_t print_host_mem();
    size_t print_device_mem();
#endif

public:
    /**
     * The constructor of StaticBTESolver
     * @param mesh      a pointer to BTEMesh
     * @param bcs       a pointer to BTEBoundaryCondition
     * @param bands     a pointer to BTEBand
     *
     * All three pointer should not be null.
     */
    StaticBTESolver(BTEMesh* mesh, BTEBoundaryCondition* bcs, BTEBand* bands);
    /**
     * Set solver hyper parameter
     * @param DM            material dimension
     * @param num_theta     CADOM parameter
     * @param num_phi       CADOM parameter
     * @param T_ref         reference temperature
     */
    void setParam(int DM, int num_theta, int num_phi, double T_ref);
    /**
     * After setting the solver context by StaticBTESolver::setParam,
     * this method can be called to invoke the solver
     * @param max_iter      maximum number of iteration allowed before exit
     */
    void solve(int max_iter);
};

#endif //BTESOLVER_STATICBTESOLVER_H
