//
// Created by Yucheng Shi on 7/10/20.
//

#include <iostream>
#include <string>
#include "StaticBTESolver/StaticBTESolver.h"

using namespace std;

int main() {
    double L_x, L_y, L_z, Tref, errnum, qflux;
    int maxIter;
    bool backup;
    int ntheta, nphi;


    ifstream infile("inputdata3D.dat");
    string line;
    char newline;
    getline(infile, line);
    infile >> L_x >> newline;
    getline(infile, line);
    infile >> L_y >> newline;
    getline(infile, line);
    infile >> L_z >> newline;
    getline(infile, line);
    infile >> Tref >> newline;
    getline(infile, line);
    infile >> errnum >> newline;
    getline(infile, line);
    infile >> ntheta >> newline;
    getline(infile, line);
    infile >> nphi >> newline;
    getline(infile, line);
    infile >> maxIter >> newline;
    getline(infile, line);
    infile >> backup >> newline;
    infile.close();

    ifstream inFile1("Input-dispersion-relation-fp.dat");
    auto bands = new BTEBand(inFile1);
    ifstream inFile2("inputbc3D.dat");
    auto bcs = new BTEBoundaryCondition(inFile2);
    ifstream inFile3("mesh3D.mphtxt");
    auto mesh3D = new  BTEMesh(inFile3, L_x, L_y, L_z);

    StaticBTESolver solver(mesh3D, bcs, bands);
    solver.setParam(3, ntheta, nphi, Tref);
    solver.solve(maxIter);

    return 0;
}