//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_BTEBOUNDARYCONDITION_H
#define BTESOLVER_BTEBOUNDARYCONDITION_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "utility.h"

class BTEBoundaryCondition {
    using BoundaryConditionList = std::vector<BoundaryCondition>;
public:
    BoundaryConditionList boundaryConditions;
    BTEBoundaryCondition() = default;
    explicit BTEBoundaryCondition(std::ifstream& inFile);
    BTEBoundaryCondition(const BTEBoundaryCondition& other) = default;
    BTEBoundaryCondition& operator=(const BTEBoundaryCondition& other) = default;

    using iterator = BoundaryConditionList::iterator;
    using const_iterator = BoundaryConditionList::const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
    std::size_t size() const;

    const BoundaryCondition& operator[](std::size_t i) const;
    friend std::ostream& operator<<(std::ostream& os, const BTEBoundaryCondition& bc);
};

#endif //BTESOLVER_BTEBOUNDARYCONDITION_H