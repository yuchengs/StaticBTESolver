//
// Created by Yucheng Shi on 7/10/20.
//

#include <iostream>
#include <string>
#include "StaticBTESolver/StaticBTESolver.h"

using namespace std;

int main() {
    double L_x, L_y, Tref, errnum, qflux;
    int maxIter;
    bool backup;
    int ntheta, nphi;


    ifstream infile("inputdata2D.dat");
    string line;
    char newline;
    getline(infile, line);
    infile >> L_x >> newline;
    getline(infile, line);
    infile >> L_y >> newline;
    getline(infile, line);
    infile >> Tref >> newline;
    getline(infile, line);
    infile >> errnum >> newline;
    getline(infile, line);
    infile >> qflux >> newline;
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
    ifstream inFile2("inputbc2D.dat");
    auto bcs = new BTEBoundaryCondition(inFile2);
    ifstream inFile3("mesh2D.mphtxt");
    auto mesh2D = new  BTEMesh(inFile3, L_x, L_y);

    StaticBTESolver solver(mesh2D, bcs, bands);
    solver.setParam(2, ntheta, nphi, Tref);
    solver.solve(maxIter);

    return 0;
}