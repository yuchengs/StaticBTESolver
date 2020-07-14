//
// Created by Yucheng Shi on 7/8/20.
//

#include <iostream>
#include <fstream>
#include "BTESolver/BTEBand.h"

using namespace std;

int main() {
    ifstream inFile("Input-dispersion-relation-fp.dat");
    BTEBand bands(inFile);
    cout << bands << endl;
}