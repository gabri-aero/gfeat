#ifndef _SINEX_HPP_
#define _SINEX_HPP_

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

class Sinex {
    std::string filename;

    int l_idx(int l) { return l * l - 4; }

public:
    Sinex(std::string filename) : filename(filename) {}
    Eigen::MatrixXd get_normal_matrix(int l_max) {
        std::ifstream file(this->filename, std::ios::in);

        if (!file.is_open()) {
            std::cerr << "Error: file " << filename << " was not opened!"
                      << std::endl;
        }

        std::string line;
        line.reserve(100);
        int n = l_idx(l_max + 1);
        Eigen::MatrixXd N(n, n);
        int row, col;
        double value;
        bool normal_equations_flag = false;

        while (getline(file, line)) {
            if (normal_equations_flag) {
                // Assign to parameters
                std::stringstream ss(line);
                ss >> row >> col >> value;
                // Map to zero-indexing
                row--;
                col--;
                if (row < n) {
                    while (!ss.fail() and col < n) {
                        N.coeffRef(col, row) = value;
                        N.coeffRef(row, col++) = value;
                        ss >> value;
                    }
                } else {
                    break;
                }
                // Check if normal equation section reached
            } else {
                if (line.find("NORMAL_EQUATION_MATRIX") != std::string::npos) {
                    normal_equations_flag = true;
                    getline(file, line); // Skip header
                }
            }
        }
        return N;
    }
};

#endif // _SINEX_HPP_
