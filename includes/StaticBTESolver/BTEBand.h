//
// Created by Yucheng Shi on 7/6/20.
//

#ifndef BTESOLVER_BTEBAND_H
#define BTESOLVER_BTEBAND_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "utility.h"

using namespace staticbtesolver;

class BTEBand {
    using BandList = std::vector<Band>;
    BandList bands;
public:
    BTEBand() = default;

    /**
     * Constructor of BTEBand
     * @param inFile inFile should be an opened file stream for a file specifying band information.
     *               The file contains one line for each band. Each line contains four numbers,
     *               "GROUP_VELOCITY RELAXATION_TIME Ctot Lr"
     */
    explicit BTEBand(std::ifstream& inFile);
    BTEBand(const BTEBand& other) = default;
    BTEBand& operator=(const BTEBand& other) = default;

    using iterator = BandList::iterator;
    using const_iterator = BandList::const_iterator;

    /**
     *
     * @return the iterator of the first band
     */
    const_iterator begin() const;

    /**
     *
     * @return the iterator of the last band
     */
    const_iterator end() const;

    /**
     *
     * @return the number of bands
     */
    std::size_t size();

    const Band& operator[](std::size_t i) const;
    friend std::ostream& operator<<(std::ostream& os, const BTEBand& band);
};

#endif //BTESOLVER_BTEBAND_H
