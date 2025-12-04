#ifndef _AOD1B_HPP_
#define _AOD1B_HPP_

#include "SphericalHarmonics.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

struct DateTime {
    int year;
    int month;
    int day;
    int h;
    int m;
    int s;

    DateTime() : year(2000), month(1), day(1), h(0), m(0), s(0) {}
    DateTime(int year, int month, int day, int h = 0, int m = 0, int s = 0)
        : year(year), month(month), day(day), h(h), m(m), s(s) {}

    bool operator==(const DateTime &other) const {
        return year == other.year && month == other.month && day == other.day &&
               h == other.h && m == other.m && s == other.s;
    }
};

enum AOD1BType { ATM, OCN, GLO, OBA };

namespace std {
template <> struct hash<DateTime> {
    size_t operator()(const DateTime &dt) const {
        size_t h1 = std::hash<int>{}(dt.year);
        size_t h2 = std::hash<int>{}(dt.month);
        size_t h3 = std::hash<int>{}(dt.day);
        size_t h4 = std::hash<int>{}(dt.h);
        size_t h5 = std::hash<int>{}(dt.m);
        size_t h6 = std::hash<int>{}(dt.s);
        return (((((h1 ^ (h2 << 1)) ^ (h3 << 2)) ^ (h4 << 3)) ^ (h5 << 4)) ^
                (h6 << 5));
    }
};

template <> struct hash<std::pair<DateTime, AOD1BType>> {
    size_t operator()(const std::pair<DateTime, AOD1BType> &p) const {
        size_t h1 = std::hash<DateTime>{}(p.first);
        size_t h2 = std::hash<int>{}(static_cast<int>(p.second));
        return h1 ^ (h2 << 1); // Combine hashes
    }
};
} // namespace std

class AOD1B {
    std::unordered_map<std::pair<DateTime, AOD1BType>, SphericalHarmonics>
        sh_models;

public:
    AOD1B() {};
    AOD1B load(std::string filename) {
        std::ifstream file(filename, std::ios::in);

        if (!file.is_open()) {
            std::cerr << "Error: file was not opened!" << std::endl;
        }

        std::string line;
        line.reserve(100);
        int l_max;
        bool read_dataset = false;
        std::regex pattern(
            R"(DATA SET (\d+):\s+\d+\s+COEFFICIENTS FOR (\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2}) OF TYPE (\w+))");
        std::smatch match;
        DateTime date;
        int dataset_number = 0;
        int total_datasets;
        int l, m;
        double Clm, Slm;
        SphericalHarmonics sh;
        // Read file header
        while (getline(file, line)) {
            if (line.find("MAXIMUM DEGREE") != std::string::npos) {
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    if (token == ":") {
                        ss >> l_max;
                        sh.l_max = l_max;
                        sh.m_max = l_max;
                        sh.coefficients.resize(l_max + 1, l_max + 1);
                        break;
                    }
                }
            } else if (line.find("CONSTANT GM") != std::string::npos) {
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    if (token == ":") {
                        ss >> sh.mu;
                    }
                }
            } else if (line.find("CONSTANT A") != std::string::npos) {
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    if (token == ":") {
                        ss >> sh.ae;
                    }
                }
            } else if (line.find("NUMBER OF DATA SETS") != std::string::npos) {
                std::stringstream ss(line);
                std::string token;
                while (ss >> token) {
                    if (token == ":") {
                        ss >> total_datasets;
                    }
                }
            } else if (line.find("END OF HEADER") != std::string::npos) {
                break;
            }
        }

        while (getline(file, line)) {
            AOD1BType type;
            if (read_dataset) {
                std::stringstream ss(line);
                ss >> l >> m >> Clm >> Slm;
                sh.Clm(l, m) = Clm;
                if (m > 0)
                    sh.Slm(l, m) = Slm;
                if (l == l_max and m == l_max) { // end of data set reached
                    read_dataset = false;
                    sh_models[{date, type}] = sh;
                }
            } else if (total_datasets == dataset_number) {
                break;
            } else {
                // Read data set header
                read_dataset = true;
                std::regex_search(line, match, pattern);
                // Get dataset number
                dataset_number = std::stoi(match[1]);
                // Get datetime
                date.year = std::stoi(match[2]);
                date.month = std::stoi(match[3]);
                date.day = std::stoi(match[4]);
                date.h = std::stoi(match[5]);
                date.m = std::stoi(match[6]);
                date.s = std::stoi(match[7]);
                // Get dataset type
                if (match[8] == "atm") {
                    type = ATM;
                } else if (match[8] == "ocn") {
                    type = OCN;
                } else if (match[8] == "glo") {
                    type = GLO;
                } else if (match[8] == "oba") {
                    type = OBA;
                } else {
                    type = UNKNOWN;
                }
                // Reset coefficients data
                sh.coefficients.setZero();
            }
        }
        return *this;
    }
    SphericalHarmonics get(DateTime date, AOD1BType type) {
        return sh_models[{date, type}];
    }
};

#endif // _AOD1B_HPP_