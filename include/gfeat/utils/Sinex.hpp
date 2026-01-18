#ifndef _SINEX_HPP_
#define _SINEX_HPP_

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
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

        // Skip to start of normal matrix data
        while (getline(file, line)) {
            if (line.find("NORMAL_EQUATION_MATRIX") != std::string::npos) {
                getline(file, line);
                break;
            }
        }

        // Read normal matrix data
        while (getline(file, line)) {
            // Check if end of normal matrix reached
            if (line[0] == '-')
                break;
            // Define line data pointer
            char *p = line.data();
            // Retrieve row and column attributes
            row = std::strtol(p, &p, 10);
            col = std::strtol(p, &p, 10);
            // Map to zero-indexing
            row--;
            col--;
            // Check normal matrix size not exceeded
            if (row >= n) {
                std::cerr << "Normal matrix size exceeded" << std::endl;
                break;
            }
            // Read and assign values
            while (*p and col < n) {
                // Only fill lower triangular part
                value = std::strtod(p, &p);
                N.coeffRef(col++, row) = value;
            }
        }
        // Now fill upper triangular part
        N.triangularView<Eigen::Upper>() = N.transpose();
        return N;
    }
};

#endif // _SINEX_HPP_
