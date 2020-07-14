//
// Created by Yucheng Shi on 7/8/20.
//

#include <iostream>
#include <fstream>
#include "BTESolver/BTEBoundaryCondition.h"

using namespace std;

int main() {
    ifstream inFile("inputbc2D.dat");
    BTEBoundaryCondition bcs(inFile);
    cout << bcs << endl;
}