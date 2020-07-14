//
// Created by Yucheng Shi on 7/6/20.
//

#include "StaticBTESolver/BTEBoundaryCondition.h"

BTEBoundaryCondition::BTEBoundaryCondition(std::ifstream& inFile) {
    if (!inFile.is_open()) {
        std::cout << "DEBUG: file not open" << std::endl;
        exit(1);
    }
    int num;
    std::string line;
    inFile >> num;
    for (int i = 0; i < 3; i++)
        getline(inFile, line);
    for (int i = 0; i < num; i++) {
        BoundaryCondition bc;
        int temp;
        inFile >> temp >> bc.type >> bc.temperature;
        bc.index = - temp - 1;
        boundaryConditions.push_back(bc);
    }
}


BTEBoundaryCondition::const_iterator BTEBoundaryCondition::begin() const {
    return boundaryConditions.begin();
}

BTEBoundaryCondition::const_iterator BTEBoundaryCondition::end() const {
    return boundaryConditions.end();
}

std::size_t BTEBoundaryCondition::size() const {
    return boundaryConditions.size();
}

const BoundaryCondition& BTEBoundaryCondition::operator[](std::size_t i) const {
    return boundaryConditions[i];
};

std::ostream& operator<<(std::ostream& os, const BTEBoundaryCondition& bc) {
    for (auto& b : bc) {
        os << b << std::endl;
    }
    return os;
}