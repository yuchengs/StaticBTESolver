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
    using BoundaryConditionList = std::vector<staticbtesolver::BoundaryCondition>;
public:
    BoundaryConditionList boundaryConditions;
    BTEBoundaryCondition() = default;

    /**
     * The constructor for BTEBoundaryCondition::BTEBoundaryCondition
     * @param inFile    inFile should be an opened file stream that specifies the boundary conditions
     *                  in the following form. The file begins with one number that specifies the
     *                  number of boundaries. Next two lines is skipped. Starting from the fourth line,
     *                  each line specifies one boundary condition as follows:
     *                  "Boundary_condition_number Boundary_condition_type Boundary_condition_temperature"
     */
    explicit BTEBoundaryCondition(std::ifstream& inFile);
    BTEBoundaryCondition(const BTEBoundaryCondition& other) = default;
    BTEBoundaryCondition& operator=(const BTEBoundaryCondition& other) = default;

    using iterator = BoundaryConditionList::iterator;
    using const_iterator = BoundaryConditionList::const_iterator;

    /**
     *
     * @return iterator of the first boundary condition
     */
    const_iterator begin() const;

    /**
     *
     * @return iterator of the last boundary condition
     */
    const_iterator end() const;

    /**
     *
     * @return the number of boundary conditions
     */
    std::size_t size() const;

    const staticbtesolver::BoundaryCondition& operator[](std::size_t i) const;
    friend std::ostream& operator<<(std::ostream& os, const BTEBoundaryCondition& bc);
};

#endif //BTESOLVER_BTEBOUNDARYCONDITION_H