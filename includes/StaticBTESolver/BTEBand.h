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

class BTEBand {
    using BandList = std::vector<Band>;
    BandList bands;
public:
    BTEBand() = default;
    explicit BTEBand(std::ifstream& inFile);
    BTEBand(const BTEBand& other) = default;
    BTEBand& operator=(const BTEBand& other) = default;

    using iterator = BandList::iterator;
    using const_iterator = BandList::const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
    std::size_t size();

    const Band& operator[](std::size_t i) const;
    friend std::ostream& operator<<(std::ostream& os, const BTEBand& band);
};

#endif //BTESOLVER_BTEBAND_H
