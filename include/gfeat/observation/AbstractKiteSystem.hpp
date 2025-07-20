#ifndef _ABSTRACT_KITE_SYSTEM_
#define _ABSTRACT_KITE_SYSTEM_

#include "../gravity/Functionals.hpp"
#include "../utils/Logger.hpp"

#include <Eigen/Dense>
#include <functions>
#include <iostream>
#include <map>
#include <set>
#include <vector>

// Custom hash function for pair
struct _PairHash {
    std::size_t operator()(const std::pair<int, int> &p) const {
        // Combine the hash of both integers
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

namespace Eigen {
using DiagonalMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;
}

class AbstractKiteSystem {
protected:
    std::vector<Eigen::MatrixXd> H; // Design matrix
    std::vector<Eigen::DiagonalMatrixXd>
        Pyy; // Observation covariance matrix (no correlation)
    std::vector<Eigen::MatrixXd> Pxx; // Parameter covariance matrix
    std::vector<Eigen::DiagonalMatrixXd>
        Pcc_inv; // Parameter a-priori covariance (no cross-terms)
    std::vector<Eigen::MatrixXd> N; // Normal matrix

    Eigen::MatrixXd sigma_x; // Parameters standard deviation matrix

    int l_max;           // Maximum degree
    std::vector<int> nx; // Number of parameters (per block)
    std::vector<int> ny; // Number of observations (per block)
    int Nx;              // Total number of parameters

    // Repeating ground-track conditions
    int Nr; // Number of revolutions
    int Nd; // Nodals days

    // Frequency resolution and number of ground-track samples
    double df;
    double n_rgt = 1;

    // Index mappings
    int n_blocks; // Number of blocks
    // Get orders that are combined to form a block (m = 0, 1, 2, ..., Nr/2)
    std::vector<std::set<int>> m_blocks;
    // Clm map to block sub-idx
    std::vector<std::unordered_map<std::pair<int, int>, int, _PairHash>>
        clm_idx;
    // Slm map to block sub-idx
    std::vector<std::unordered_map<std::pair<int, int>, int, _PairHash>>
        slm_idx;
    // Clm map to global sub-idx
    std::unordered_map<std::pair<int, int>, int, _PairHash> global_clm_idx;
    // Slm map to global sub-idx
    std::unordered_map<std::pair<int, int>, int, _PairHash> global_slm_idx;
    // Global initialization flag
    bool global_initialized = false;

public:
    AbstractKiteSystem(int l_max, int Nr, int Nd)
        : l_max(l_max), Nr(Nr), Nd(Nd) {
        auto t1 = std::chrono::high_resolution_clock::now();
        // Block numbers solely defined by number of revolutions of repeating
        // orbit
        n_blocks = std::min(Nr / 2 + 1, l_max + 1);
        // Allocate everything for computed number of blocks
        m_blocks.resize(n_blocks);
        clm_idx.resize(n_blocks);
        slm_idx.resize(n_blocks);
        H.resize(n_blocks);
        Pyy.resize(n_blocks);
        Pxx.resize(n_blocks);
        Pcc_inv.resize(n_blocks);
        N.resize(n_blocks);
        nx.resize(n_blocks);
        ny.resize(this->n_blocks);
        // Define orders that form a common block
        int n, m;
        // Loop over blocks
        for (int m0 = 0; m0 < n_blocks; m0++) {
            // Cache references
            auto &m_block = m_blocks[m0];
            auto &clm = clm_idx[m0];
            auto &slm = slm_idx[m0];
            // Assign orders that form every block
            n = 0;
            while (true) {
                m = n * Nr + m0;
                if (m > this->l_max)
                    break;
                m_block.insert(m);
                m = (n + 1) * Nr - m0;
                if (m > this->l_max)
                    break;
                m_block.insert(m);
                n++;
            }
            // Map Clm, Slm to sub-block index
            int block_x_idx = 0;
            for (int m : m_block) {         // Loop over orders of that block
                int l_min = std::max(m, 2); // Pre-compute lower degree limit
                // First Clm terms
                for (int l = l_min; l <= l_max; l++) {
                    clm[std::make_pair(l, m)] = block_x_idx++;
                }
                // Then Slm terms
                if (m != 0) {
                    for (int l = l_min; l <= this->l_max; l++) {
                        slm[std::make_pair(l, m)] = block_x_idx++;
                    }
                }
            }
            // Assign number of block parameters
            nx[m0] = block_x_idx;
            // Allocate parameter matrices
            Pxx[m0].resize(nx[m0], nx[m0]);
            Pcc_inv[m0].resize(nx[m0]);
            N[m0].resize(nx[m0], nx[m0]);
            // Set to zero
            Pxx[m0].setZero();
            Pcc_inv[m0].setZero();
            N[m0].setZero();
        }
        // Allocate sigma_x
        sigma_x.resize(l_max + 1, l_max + 1);
        sigma_x.setZero();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        Logger::instance() << "Block system - parameter initialization: "
                           << duration.count() * 1e-3 << " s\n";
    }
    void global_initialization() {
        // Define typical global indexing for parameters
        int i = 0;
        for (int m = 0; m <= l_max; m++) {
            for (int l = std::max(m, 2); l <= l_max; l++) {
                global_clm_idx[std::make_pair(l, m)] = i++;
            }
            for (int l = std::max(m, 2); l <= l_max; l++) {
                if (m > 0)
                    global_slm_idx[std::make_pair(l, m)] = i++;
            }
        }
        // Compute full dimensions
        this->Nx = 0;
        for (int nx : this->nx) {
            this->Nx += nx;
        }
        // Global initialization flag
        global_initialized = true;
    }
    void set_kaula_regularization() {
        double sigma_lm;
        for (int m0 = 0; m0 < n_blocks; m0++) {
            for (int m : m_blocks.at(m0)) {
                for (int l = std::max(m, 2); l <= l_max; l++) {
                    sigma_lm = 1e-5 / pow(l, 2); // Kaula's rule of thumb
                    int clm_i = clm_idx.at(m0)[{l, m}];
                    Pcc_inv.at(m0).diagonal()(clm_i) = 1 / pow(sigma_lm, 2);
                    if (m > 0) {
                        int slm_i = slm_idx.at(m0)[{l, m}];
                        Pcc_inv.at(m0).diagonal()(slm_i) = 1 / pow(sigma_lm, 2);
                    }
                }
            }
        }
    }

    void block_solve() {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int m0 = 0; m0 < n_blocks; m0++) {
            // Compute weighted normal matrix
            N.at(m0) = H.at(m0).transpose() *
                       (1 / n_rgt * Pyy.at(m0))
                           .diagonal()
                           .cwiseInverse()
                           .asDiagonal() *
                       H.at(m0);
            // Compute parameter covariance
            // Cholesky factorization employed for efficient inversion
            Pxx.at(m0) = (N.at(m0) + Pcc_inv.at(m0).toDenseMatrix())
                             .llt()
                             .solve(Eigen::MatrixXd::Identity(N.at(m0).rows(),
                                                              N.at(m0).cols()));
            // Populate sigma_x matrix
            for (auto &[lm, i] : clm_idx.at(m0)) {
                auto [l, m] = lm;
                sigma_x(l, m) = sqrt(Pxx.at(m0)(i, i));
            }
            for (auto &[lm, i] : slm_idx.at(m0)) {
                auto [l, m] = lm;
                if (m > 0)
                    sigma_x(m - 1, l) = sqrt(Pxx.at(m0)(i, i));
            }
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        Logger::instance() << "Parameter covariance solved: "
                           << duration.count() * 1e-3 << " s\n";
    }
    auto &get_Pyy_blocks() { return Pyy; }
    auto &get_Pxx_blocks() { return Pxx; }
    auto &get_N_blocks() { return N; }
    auto &get_H_blocks() { return H; }
    auto &get_sigma_x() { return sigma_x; }
    auto get_nx() const { return nx; }
    auto get_ny() const { return ny; }
    virtual Eigen::MatrixXd get_H() = 0;
    virtual Eigen::MatrixXd get_Pyy() = 0;
    Eigen::MatrixXd get_N() {
        Eigen::MatrixXd H = get_H();
        Eigen::MatrixXd Pyy = get_Pyy();
        Eigen::MatrixXd N =
            H.transpose() * Pyy.diagonal().cwiseInverse().asDiagonal() * H;
        return N;
    }
    double get_sigma_Clm(int l, int m) const { return sigma_x(l, m); }
    double get_sigma_Slm(int l, int m) const {
        return m == 0 ? 0 : sigma_x(m - 1, l);
    }
    auto degree_variance() const {
        // Allocate degree variance vector
        Eigen::VectorXd sigma2_l(l_max + 1);
        sigma2_l.setZero();
        for (int l = 2; l <= l_max; l++) {
            // Clm terms
            for (int m = 1; m <= l; m++) {
                sigma2_l[l] += pow(this->get_sigma_Clm(l, m), 2);
            }
            // Slm terms
            for (int m = 1; m <= l; m++) {
                sigma2_l[l] += pow(this->get_sigma_Slm(l, m), 2);
            }
        }
        return sigma2_l;
    }

    Eigen::VectorXd rms_per_coefficient_per_degree() const {
        Eigen::VectorXd result(l_max + 1);
        auto deg_var = this->degree_variance();
        for (int l = 2; l <= this->l_max; l++) {
            result(l) = std::sqrt(deg_var.coeff(l) / (2 * l + 1));
        }
        return result;
    }

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
    synthesis(int n_lon, int n_lat, const BaseFunctional &functional) const {
        auto t1 = std::chrono::high_resolution_clock::now();
        // Create longitude and latitude vectors (at center of discrete
        // intervals)
        double dlon = 2 * M_PI / n_lon;
        double dlat = M_PI / n_lat;
        Eigen::VectorXd lon_vec =
            Eigen::VectorXd::LinSpaced(n_lon, -n_lon / 2 + 1, n_lon / 2) * dlon;
        Eigen::VectorXd colat_vec =
            Eigen::VectorXd::LinSpaced(n_lat, 0.5, n_lat - 0.5) * dlat;
        Eigen::VectorXd lat_vec = (colat_vec.array() - M_PI_2).matrix();
        // Apply broadcasting for generating mesh
        Eigen::MatrixXd lon_mesh =
            (lon_vec * 180 / M_PI).transpose().replicate(n_lat, 1);
        Eigen::MatrixXd lat_mesh = (lat_vec * 180 / M_PI).replicate(1, n_lon);
        // Define functional std deviation
        Eigen::MatrixXd y(n_lat, n_lon);
        y.setZero();
        // Define cos_mlam, sin_mlam
        Eigen::MatrixXd cos_mlam(this->l_max + 1, n_lon),
            sin_mlam(this->l_max + 1, n_lon);
        for (int lon_idx = 0; lon_idx < n_lon; lon_idx++) {
            for (int m = 0; m <= this->l_max; m++) {
                cos_mlam(m, lon_idx) = cos(m * lon_vec(lon_idx));
                sin_mlam(m, lon_idx) = sin(m * lon_vec(lon_idx));
            }
        }
        // Define f constant
        Eigen::VectorXd f(this->l_max + 1);
        for (int l = 0; l <= this->l_max; l++) {
            f(l) = functional.degree_dependent_factor(l);
        }
        f = f * functional.common_factor;
        // Define Am, Bm
        Eigen::VectorXd A(this->l_max + 1), B(this->l_max + 1);
        // Define latitude loop
        for (int lat_idx = 0; lat_idx < n_lat; lat_idx++) {
            Logger::instance() << "Propagating covariance: " << lat_idx + 1
                               << "/" << n_lat << "\r" << std::flush;
            // Compute Legendre polynomials
            Plm plm(l_max, colat_vec(lat_idx));
            // Loop over block covariance matrices
            for (int m0 = 0; m0 < n_blocks; m0++) {
                // Get block references
                auto &Pxx = this->Pxx[m0];
                auto clm = this->clm_idx[m0];
                auto slm = this->slm_idx[m0];
                int nx_block = this->nx[m0];
                // Pre-allocate vector of partials and W matrix
                Eigen::MatrixXd W(nx_block, n_lon);
                Eigen::VectorXd v(nx_block);
                v.setZero();
                W.setZero();
                // Loop over block-covariance rows
                for (int row_idx = 0; row_idx < nx_block; row_idx++) {
                    // Fill A, B (latitude dependent terms)
                    A.setZero();
                    B.setZero();
                    // Assign Am terms first
                    for (auto &[lm, col_idx] : clm) {
                        auto [l, m] = lm;
                        A(m) += f(l) * plm.get_Plm_bar(l, m) *
                                Pxx(row_idx, col_idx);
                    }
                    // Then Bm terms
                    for (auto &[lm, col_idx] : slm) {
                        auto [l, m] = lm;
                        B(m) += f(l) * plm.get_Plm_bar(l, m) *
                                Pxx(row_idx, col_idx);
                    }
                    // Compute Wi matrix (way more faster than FFT - Am, Bm full
                    // of zeros)
                    for (int lon_idx = 0; lon_idx < n_lon; lon_idx++) {
                        for (int m : m_blocks[m0]) {
                            W(row_idx, lon_idx) += A(m) * cos_mlam(m, lon_idx) +
                                                   B(m) * sin_mlam(m, lon_idx);
                        }
                    }
                }
                // Loop over longitudes
                for (int lon_idx = 0; lon_idx < n_lon; lon_idx++) {
                    // Fill vector of partials
                    // Assign vector of block derivatives - Clm terms first
                    for (auto &[lm, i] : clm) {
                        auto [l, m] = lm;
                        v(i) =
                            f(l) * plm.get_Plm_bar(l, m) * cos_mlam(m, lon_idx);
                    }
                    // Then Slm terms
                    for (auto &[lm, i] : slm) {
                        auto [l, m] = lm;
                        v(i) =
                            f(l) * plm.get_Plm_bar(l, m) * sin_mlam(m, lon_idx);
                    }
                    // Compute std deviation
                    y(lat_idx, lon_idx) += v.transpose() * W.col(lon_idx);
                }
            }
        }
        std::cout << "\n";
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        Logger::instance() << "Covariance propagated: "
                           << duration.count() * 1e-3 << "\n";
        return {lon_mesh, lat_mesh, y.array().sqrt()};
    }
    double synthesis_average(const BaseFunctional &functional) const {
        double mean_y_error = 0.0;
        double sigma2 = 0;
        for (int l = 2; l <= this->l_max; l++) {
            double sigma2_m = 0.0;
            for (int m = 0; m <= l; m++) {
                sigma2_m += m == 0 ? pow(this->get_sigma_Clm(l, m), 2)
                                   : pow(this->get_sigma_Clm(l, m), 2) +
                                         pow(this->get_sigma_Slm(l, m), 2);
            }
            sigma2 += pow(functional.degree_dependent_factor(l), 2) * sigma2_m;
        }
        mean_y_error = functional.common_factor * sqrt(sigma2);
        return mean_y_error;
    }

    void set_solution_time_window(double days) {
        double time_window = days * 24 * 3600;
        this->n_rgt = time_window * this->df;
        std::cout << this->n_rgt << std::endl;
    }
};

#endif //_ABSTRACT_KITE_SYSTEM_