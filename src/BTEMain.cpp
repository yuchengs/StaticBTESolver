//
// Created by Yucheng Shi on 7/13/20.
//

#include <iostream>
#include <string>
#include <getopt.h>
#include "StaticBTESolver/BTEGeometry.h"
#include "StaticBTESolver/StaticBTESolver.h"
#ifdef USE_GPU
#include <mpi.h>
#endif

using namespace std;

/*
 * Command Line Utility to use the solver
 * Usage: -g    follow by a geo file, 3 double representing Lx, Ly, Lz
 *        -m    follow by a .dat file specifying bands information
 *        -b    follow by a .dat file specifying each boundary condition,
 *              if this argument is not provided, you have to provide boundary conditions
 *              interactively in the command line
 */

int main (int argc, char **argv) {
    MPI_Init(nullptr, nullptr);
    opterr = true;
    static struct option longopts[] = {
            { "geometry", required_argument, nullptr, 'g' },
            { "material", required_argument, nullptr, 'm' },
            { "boundary", optional_argument, nullptr, 'b' },
            { "DM", required_argument, nullptr, 'd' },
            { "ntheta", required_argument, nullptr, 't' },
            { "nphi", required_argument, nullptr, 'p' },
            { "WFACTOR", optional_argument, nullptr, 'w' },
            { "T_ref", optional_argument, nullptr, 'T' },
            { "maxIter", optional_argument, nullptr, 'I' },
            { "L_x", optional_argument, nullptr, 'x' },
            { "L_y", optional_argument, nullptr, 'y' },
            { "L_z", optional_argument, nullptr, 'z' },


            { nullptr, 0, nullptr, '\0' }
    };

    string geofileName, bandfileName, bfileName;
    double L_x = 0, L_y = 0, L_z = 0;

    int DM = 3, ntheta = 4, nphi = 4;
    double T_ref = 300, WFACTOR = 1;
    int maxIter = 10000;

    int idx = 0;
    int c;
    bool bPresent = false;

    while ((c = getopt_long(argc, argv, "g:m:b:d:t:p:w:T:I:x:y:z:", longopts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                string str(optarg);
                stringstream ss(str);
                ss >> geofileName;
                break;
            }
            case 'm': {
                bandfileName = optarg;
                break;
            }
            case 'b': {
                bPresent = true;
                bfileName = optarg;
                break;
            }
            case 'd': {
                DM = stoi(optarg);
                break;
            }
            case 't': {
                ntheta = stoi(optarg);
                break;
            }
            case 'p': {
                nphi = stoi(optarg);
                break;
            }
            case 'w': {
                std::string str(optarg);
                if (str.back() >= '0' && str.back() <= '9') {
                    WFACTOR = stod(str);
                }
                else {
                    WFACTOR = stod(str) * PI;
                }
                break;
            }
            case 'T': {
                T_ref = stod(optarg);
                break;
            }
            case 'I': {
                maxIter = stoi(optarg);
                break;
            }
            case 'x': {
                L_x = stod(optarg);
                break;
            }
            case 'y': {
                L_y = stod(optarg);
                break;
            }
            case 'z': {
                L_z = stod(optarg);
                break;
            }
            default: {
                cerr << "Unknown option " << c << endl;
                return 1;
            }
        }
    }

    BTEMesh* mesh;
    BTEGeometry* geo;
    string ext;
    int N_cell = 0;
    auto p = geofileName.find_last_of('.');
    if (p == string::npos) {
        N_cell = stoi(geofileName);
        mesh = new BTEMesh(N_cell, L_x);
    }
    else {
        ext = geofileName.substr(p + 1);
        if (ext == "geo") {
            geo = new BTEGeometry(geofileName, L_x, L_y, L_z);
            mesh = geo->export_mesh();
        } else if (ext == "mphtxt") {
            ifstream geofile(geofileName);
            mesh = new BTEMesh(geofile, L_x, L_y, L_z);
            geofile.close();
        }
        else {
            cout << "file format not supported" << endl;
            exit(1);
        }
    }

    ifstream bandFile(bandfileName);
    auto bands = new BTEBand(bandFile);
    bandFile.close();
    BTEBoundaryCondition* bcs;
    if (bPresent) {
        ifstream bcFile(bfileName);
        bcs = new BTEBoundaryCondition(bcFile);
        bcFile.close();
    }
    else {
        if (p == string::npos || ext != "geo") {
            cerr << "Must specify boundary conditions" << endl;
            exit(1);
        }
        cout << "You have not specified boundary conditions for ";
        cout << geo->boundary_indices.size() << " boundaries." << endl;
        bcs = new BTEBoundaryCondition();
        bcs->boundaryConditions.reserve(geo->boundary_indices.size());
        for (int boundary_index : geo->boundary_indices) {
            cout << "Boundary #" << boundary_index << ":" << endl;
            BoundaryCondition bc(- boundary_index - 1);
            cout << "Type: ";
            cin >> bc.type;
            cout << "Temperature: ";
            cin >> bc.temperature;
            cout << endl;
            bcs->boundaryConditions.push_back(bc);
        }
    }

    StaticBTESolver solver(mesh, bcs, bands);
    solver.setParam(DM, ntheta, nphi, WFACTOR, T_ref);
    solver.solve(maxIter);

    delete mesh;
    delete bcs;
    delete bands;
    MPI::Finalize();
}