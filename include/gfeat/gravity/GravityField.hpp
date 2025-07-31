#ifndef _GRAVITY_FIELD_HPP_
#define _GRAVITY_FIELD_HPP_

#include "../utils/FileReader.hpp"
#include "SphericalHarmonics.hpp"

class GravityField : public SphericalHarmonics {
public:
    using SphericalHarmonics::SphericalHarmonics;

    GravityField load(std::string filename) {
        // Read data
        auto data = FileReader::read(filename, ' ');
        int l, m;
        // Assign data
        for (auto row : data) {
            if (row[0] == "earth_gravity_constant") {
                // Handle Fortran D for exponents
                std::string value_str = row[1];
                size_t pos = value_str.find_first_of("Dd");
                if (pos != std::string::npos) {
                    value_str.replace(pos, 1, "E");
                }
                this->mu = std::stod(value_str);
            } else if (row[0] == "radius") {
                // Handle Fortran D for exponents
                std::string value_str = row[1];
                size_t pos = value_str.find_first_of("Dd");
                if (pos != std::string::npos) {
                    value_str.replace(pos, 1, "E");
                }
                this->ae = std::stod(value_str);
            } else if (row[0] == "gfc" || row[0] == "gfct") {
                l = std::stoi(row[1]);
                m = std::stoi(row[2]);
                if (l <= this->l_max and m <= this->l_max) {
                    this->Clm(l, m) = std::stod(row[3]);
                    if (m > 0)
                        this->Slm(l, m) = std::stod(row[4]);
                }
            }
        }
        return *this;
    }
};

#endif // _GRAVITY_FIELD_HPP_