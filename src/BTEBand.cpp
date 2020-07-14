//
// Created by Yucheng Shi on 7/6/20.
//

#include <sstream>
#include "StaticBTESolver/BTEBand.h"

BTEBand::BTEBand(std::ifstream& inFile) {
    if (!inFile.is_open()) {
        std::cout << "DEBUG: file not open" << std::endl;
        exit(1);
    }
    std::string line;
    while (getline(inFile, line)) {
        if (line == "\r") {
            break;
        }
        std::stringstream ss(line);
        Band b;
        ss >> b.group_velocity
           >> b.relaxation_time
           >> b.Ctot
           >> b.Lr;
        bands.push_back(b);
    }
}

BTEBand::const_iterator BTEBand::begin() const {
    return bands.begin();
}

BTEBand::const_iterator BTEBand::end() const {
    return bands.end();
}

std::size_t BTEBand::size() {
    return bands.size();
}

const Band& BTEBand::operator[](std::size_t i) const {
    return bands[i];
}

std::ostream& operator<<(std::ostream& os, const BTEBand& band) {
    for (auto& b : band) {
        os << b << std::endl;
    }
    return os;
}