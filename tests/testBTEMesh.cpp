//
// Created by Yucheng Shi on 7/8/20.
//

#include <iostream>
#include <fstream>
#include "StaticBTESolver/BTEMesh.h"

using namespace std;

int main() {
    ifstream inFile2D("mesh2D.mphtxt");
    BTEMesh mesh2D(inFile2D, 1, 1);

    ifstream inFile3D("mesh3D.mphtxt");
    BTEMesh mesh3D(inFile3D, 1, 1, 1);
}