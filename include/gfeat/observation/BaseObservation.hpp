#ifndef _BASE_OBSERVATION_HPP_
#define _BASE_OBSERVATION_HPP_

#include "../gravity/GravityField.hpp"
#include "../utils/Logger.hpp"
#include "../utils/RGT.hpp"
#include "AbstractKiteSystem.hpp"

#include <functional>
#include <functions>

class BaseObservation : public AbstractKiteSystem {
protected:
    std::vector<Eigen::MatrixXd> H; // Design matrix
    std::vector<Eigen::DiagonalMatrixXd>
        Pyy; // Observation covariance matrix (no correlation)

    std::vector<int> ny; // Number of observations (per block)
    int Ny;              // Total number of observations
    // Index mappings
    // (k,m) to observation block index mapping
    std::vector<std::unordered_map<std::pair<int, int>, int, _PairHash>>
        Akm_idx;
    std::vector<std::unordered_map<std::pair<int, int>, int, _PairHash>>
        Bkm_idx;
    // Global observation indexing (Akm, Bkm)
    std::unordered_map<std::pair<int, int>, int, _PairHash> global_Akm_idx;
    std::unordered_map<std::pair<int, int>, int, _PairHash> global_Bkm_idx;
    // (Integer) frequencies per block
    std::vector<std::set<int>> gamma_km; // k*Nr - m*Nd

    // Configuration parameters
    double r;                          // Orbital radius
    double I;                          // Inclination
    std::function<double(double)> asd; // Amplitude Spectral Density
    double n;                          // Mean motion
    double we_dot;                     // Earth frequency
    double wo_dot;                     // Orbital frequency
    double we_0; // Initial Earth angle (initial longitude)
    double wo_0; // Initial argument of latitude
    Flmp flmp;   // Inclination functions

    // Shift in initial phase angle to handle general
    // purpose observation
    // 0    : ( Clm,  Slm)
    // π/2  : (-Slm, Clm)
    // π    : (-Clm, -Slm)
    // -π/2 : ( Slm, -Clm)
    double delta_psi_km_0 = 0;

    // Function that defines partials
    virtual double H_(int l, int m, int k) const = 0;

public:
    // Function that defines starting cut-off degree for summation
    virtual int l_min(int m, int k) const {
        int delta = (k - std::max({std::abs(k), m, 2})) % 2 == 0 ? 0 : 1;
        return std::max({std::abs(k), m, 2}) + delta;
    }

    BaseObservation(int l_max, int Nr, int Nd, double I, double we_0 = 0.0,
                    double wo_0 = 0.0)
        : AbstractKiteSystem(l_max, Nr, Nd),
          r(sma_from_rgt(Nr, Nd, 7000.0e3, I)), I(I), we_0(we_0), wo_0(wo_0) {
        // Compute mean motion
        n = sqrt(mu / pow(r, 3));
        this->df = n / (Nr * 2 * M_PI);
        // Compute wo_dot, we_dot
        we_dot = raan_dot(r, 0, I) - theta_dot;
        wo_dot = aop_dot(r, 0, I) + ma_dot(r, 0, I);
    }

    void initialize_block_observation_system() {
        auto t1 = std::chrono::high_resolution_clock::now();
        // Allocate everything observation related for computed number of blocks
        H.resize(this->n_blocks);
        Pyy.resize(this->n_blocks);
        ny.resize(this->n_blocks);
        Akm_idx.resize(this->n_blocks);
        Bkm_idx.resize(this->n_blocks);
        gamma_km.resize(this->n_blocks);
        // Loop over blocks
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            // Cache references
            auto &akm = Akm_idx[m0];
            auto &bkm = Bkm_idx[m0];
            auto &gamma = gamma_km[m0];
            auto &m_block = this->m_blocks[m0];
            // Compute sorted set of (integer) frequencies of the block
            int kNr_mNd;
            for (int k = -this->l_max; k <= this->l_max; k++) {
                for (int m : m_block) { // Loop over order of that block
                    kNr_mNd = std::abs(k * this->Nr - m * this->Nd);
                    if (kNr_mNd != 0 and kNr_mNd != this->Nr and
                        l_min(m, k) <=
                            this->l_max) { // Check non-resonant solution and
                                           // limiting frequencies
                        gamma.insert(kNr_mNd);
                    }
                }
            }
            // Assign number of block observations
            this->ny[m0] =
                2 * gamma.size(); // 2 observations per frequency [Akm, Bkm]
            // Compute (k,m) observation index mapping to block observation
            // index
            int yi;
            for (int m : m_block) { // Loop over order of that block
                for (int k = -this->l_max; k <= this->l_max; k++) {
                    kNr_mNd = std::abs(k * this->Nr - m * this->Nd);
                    auto it = gamma.find(kNr_mNd);
                    auto km = std::make_pair(k, m);
                    if (it != gamma.end()) { // Avoid resonant solutions
                        yi = std::distance(gamma.begin(), it);
                        akm[km] = 2 * yi;
                        bkm[km] = 2 * yi + 1;
                    } else {
                        akm[km] = -1; // Invalid case: global idx -1
                        bkm[km] = -1;
                    }
                }
            }
            // Allocate remaining matrices
            this->H[m0].resize(this->ny[m0], this->nx[m0]);
            this->Pyy[m0].resize(this->ny[m0]);
            // Set to zero
            this->H[m0].setZero();
            this->Pyy[m0].setZero();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        logger << "Block system - observation initialization: "
               << duration.count() * 1e-3 << " s\n";
        this->populate_design_matrix();
    }

    void global_initialization() {
        // Delegate function
        AbstractKiteSystem::global_initialization();
        // Define global set of observation frequencies
        std::set<int> gamma_km_full;
        for (const auto &s : gamma_km) {
            gamma_km_full.insert(s.begin(), s.end());
        }
        // Define global observation indexing
        int yi;
        for (int k = -this->l_max; k <= this->l_max; k++) {
            for (int m = 0; m <= this->l_max;
                 m++) { // Loop over order of that block
                int kNr_mNd = std::abs(k * this->Nr - m * this->Nd);
                auto it = gamma_km_full.find(kNr_mNd);
                auto km = std::make_pair(k, m);
                if (it != gamma_km_full.end()) { // Avoid resonant solutions
                    yi = std::distance(gamma_km_full.begin(), it);
                    global_Akm_idx[km] = 2 * yi;
                    global_Bkm_idx[km] = 2 * yi + 1;
                } else {
                    global_Akm_idx[km] = -1; // Resonant case: global idx -1
                    global_Bkm_idx[km] = -1;
                }
            }
        }
        // Compute full dimensions
        this->Ny = 0;
        for (int ny : this->ny) {
            this->Ny += ny;
        }
    }

    virtual void compute_flmp() {
        auto t1 = std::chrono::high_resolution_clock::now();
        flmp = Flmp(this->l_max, this->I, false);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        logger << "Inclination functions computed: " << duration.count() * 1e-3
               << " s\n";
    }

    auto &get_Pyy_blocks() { return Pyy; }
    auto &get_H_blocks() { return H; }
    auto get_ny() const { return ny; }

    Eigen::MatrixXd get_H() {
        if (!this->global_initialized)
            this->global_initialization();
        // Pre-allocate
        Eigen::MatrixXd H_full(this->Ny, this->Nx);
        H_full.setZero();
        // Assign from block matrices
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            for (int m : this->m_blocks.at(m0)) {
                for (int k = -this->l_max; k <= this->l_max;
                     k++) { // Loop over k-index
                    // Compute observation index
                    int block_Akm_i = Akm_idx.at(m0)[{k, m}];
                    int block_Bkm_i = Bkm_idx.at(m0)[{k, m}];
                    int global_Akm_i = global_Akm_idx[{k, m}];
                    int global_Bkm_i = global_Bkm_idx[{k, m}];
                    // Avoid resonant solutions (block observation idx set to
                    // -1)
                    if (block_Akm_i != -1) {
                        // Compute lower degree limit
                        int delta = (k - std::max({std::abs(k), m, 2})) % 2 == 0
                                        ? 0
                                        : 1;
                        int l_min = std::max({std::abs(k), m, 2}) + delta;
                        // Lump terms over the degree
                        for (int l = l_min; l <= this->l_max; l += 2) {
                            // Compute Clm, Slm indexing
                            int block_clm_i = this->clm_idx.at(m0)[{l, m}];
                            int block_slm_i = this->slm_idx.at(m0)[{l, m}];
                            int global_clm_i = this->global_clm_idx[{l, m}];
                            int global_slm_i = this->global_slm_idx[{l, m}];
                            H_full(global_Akm_i, global_clm_i) =
                                this->H.at(m0)(block_Akm_i, block_clm_i);
                            H_full(global_Bkm_i, global_clm_i) =
                                this->H.at(m0)(block_Bkm_i, block_clm_i);
                            if (m > 0) {
                                H_full(global_Bkm_i, global_slm_i) =
                                    this->H.at(m0)(block_Bkm_i, block_slm_i);
                                H_full(global_Akm_i, global_slm_i) =
                                    this->H.at(m0)(block_Akm_i, block_slm_i);
                            }
                        }
                    }
                }
            }
        }
        // Return full design matrix
        return H_full;
    }

    Eigen::MatrixXd get_Pyy() {
        if (!this->global_initialized)
            this->global_initialization();
        Eigen::DiagonalMatrixXd Pyy_full(Ny);
        Pyy_full.setZero();
        // Assign from block matrices
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            for (int m : this->m_blocks.at(m0)) {
                for (int k = -this->l_max; k <= this->l_max;
                     k++) { // Loop over k-index
                    // Compute observation index
                    int block_Akm_i = Akm_idx.at(m0)[{k, m}];
                    int block_Bkm_i = Bkm_idx.at(m0)[{k, m}];
                    int global_Akm_i = global_Akm_idx[{k, m}];
                    int global_Bkm_i = global_Bkm_idx[{k, m}];
                    // Assign values
                    if (global_Akm_i != -1)
                        Pyy_full.diagonal()(global_Akm_i) =
                            this->Pyy.at(m0).diagonal()(block_Akm_i);
                    if (global_Akm_i != -1)
                        Pyy_full.diagonal()(global_Bkm_i) =
                            this->Pyy.at(m0).diagonal()(block_Bkm_i);
                }
            }
        }
        return Pyy_full.toDenseMatrix();
    }

    void set_observation_error(std::function<double(double)> asd) {
        this->asd = asd;
        // Populate observation covariance
        double f, psd;
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            auto &Akm_idx = this->Akm_idx[m0];
            auto &Bkm_idx = this->Bkm_idx[m0];
            auto &Pyy = this->Pyy[m0];
            for (auto &[km, i] : Akm_idx) {
                if (i != -1) {
                    auto [k, m] = km;
                    // Compute frequency, PSD
                    f = abs(k * wo_dot + m * we_dot) / (2 * M_PI);
                    psd = pow(asd(f), 2) * this->df;
                    // Allocate to Akm
                    Pyy.diagonal()(i) = psd;
                }
            }
            for (auto &[km, i] : Bkm_idx) {
                if (i != -1) {
                    auto [k, m] = km;
                    // Compute frequency, PSD
                    f = abs(k * wo_dot + m * we_dot) / (2 * M_PI);
                    psd = pow(asd(f), 2) * this->df;
                    // Allocate to Bkm
                    Pyy.diagonal()(i) = psd;
                }
            }
        }
    }

    void populate_design_matrix() {
        // Compute inclination functions if required
        if (this->flmp.get_l_max() == 0)
            this->compute_flmp();
        auto t1 = std::chrono::high_resolution_clock::now();
        // Allocate variables for loop
        double Hlmk, beta_km;
        int Akm_i, Bkm_i, clm_i, slm_i;
        for (int m0 = 0; m0 < this->n_blocks;
             m0++) { // Loop over order (block-based)
                     // Assign references for block
            auto &H = this->H.at(m0);
            auto &Akm_idx = this->Akm_idx.at(m0);
            auto &Bkm_idx = this->Bkm_idx.at(m0);
            auto &clm_idx = this->clm_idx.at(m0);
            auto &slm_idx = this->slm_idx.at(m0);
            // Loop over orders in block
            for (int m : this->m_blocks.at(m0)) {
                for (int k = -this->l_max; k <= this->l_max;
                     k++) { // Loop over k-index
                    // Compute initial phase angle
                    double psi_km_0 = delta_psi_km_0 + k * wo_0 + m * we_0;
                    // Compute normalised frequency [cpr]
                    beta_km = k - m * static_cast<double>(this->Nd) / this->Nr;
                    // Compute observation index
                    Akm_i = Akm_idx[{k, m}];
                    Bkm_i = Bkm_idx[{k, m}];
                    double cos_psi_km_0 = cos(psi_km_0);
                    double sin_psi_km_0 = sin(psi_km_0);
                    // Avoid resonant solutions (block observation idx set to
                    // -1)
                    if (Akm_i != -1) {
                        for (int l = this->l_min(m, k); l <= this->l_max;
                             l += 2) {
                            // Lump terms over the degree
                            Hlmk = this->H_(
                                l, m, k); // compute design matrix common term
                            // Compute Clm, Slm indexing
                            clm_i = clm_idx[{l, m}];
                            if (m > 0)
                                slm_i = slm_idx[{l, m}];
                            // Compute sign of frequency
                            int s = beta_km > 0 ? 1 : -1;
                            // Assign Hlmk to global design matrix
                            // Assign dSlm only for m>0
                            if ((l - m) % 2 == 0) { // l-m: even
                                H(Akm_i, clm_i) +=
                                    Hlmk * cos_psi_km_0; // dClm * cos -> Akm
                                H(Bkm_i, clm_i) +=
                                    -s * Hlmk *
                                    sin_psi_km_0; // - dClm * sin -> Bkm
                                if (m > 0) {
                                    H(Akm_i, slm_i) +=
                                        Hlmk *
                                        sin_psi_km_0; // dSlm * sin -> Akm
                                    H(Bkm_i, slm_i) +=
                                        s * Hlmk *
                                        cos_psi_km_0; // dSlm * cos -> Bkm
                                }
                            } else { // l-m: odd
                                if (m > 0) {
                                    H(Akm_i, slm_i) +=
                                        -Hlmk *
                                        cos_psi_km_0; // -dSlm * cos -> Akm
                                    H(Bkm_i, slm_i) +=
                                        s * Hlmk *
                                        sin_psi_km_0; // dSlm * sin -> Bkm
                                }
                                H(Akm_i, clm_i) +=
                                    Hlmk * sin_psi_km_0; // dClm * sin -> Akm
                                H(Bkm_i, clm_i) +=
                                    s * Hlmk *
                                    cos_psi_km_0; // dClm * cos -> Bkm
                            }
                        }
                    }
                }
            }
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        logger << "Design matrix populated: " << duration.count() * 1e-3
               << " s\n";
    }
    double get_radius() const { return r; }
    double get_I() const { return I; }
    double get_wo_0() const { return wo_0; }
    double get_we_0() const { return we_0; }

    auto simulate_observations(const GravityField &gravity_field) {
        // Define block parameter and observation vectors
        std::vector<Eigen::VectorXd> x(this->n_blocks);
        std::vector<Eigen::VectorXd> y(this->n_blocks);
        // Allocate output
        std::unordered_map<int, std::pair<double, double>> output;
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            // Define references
            auto &H = this->H.at(m0);
            auto &xm = x.at(m0);
            auto &ym = y.at(m0);
            auto &clm_idx = this->clm_idx.at(m0);
            auto &slm_idx = this->slm_idx.at(m0);
            auto &akm_idx = this->Akm_idx.at(m0);
            auto &bkm_idx = this->Bkm_idx.at(m0);
            // Resize parameter vector
            xm.resize(this->nx.at(m0));
            // Loop over slm idxs
            for (auto &[lm, i] : slm_idx) {
                auto [l, m] = lm;
                xm(i) = m > 0 ? gravity_field.Slm(l, m) : 0;
            }
            // Loop over clm idxs
            for (auto &[lm, i] : clm_idx) {
                auto [l, m] = lm;
                xm(i) = gravity_field.Clm(l, m);
            }
            // Generate observations for current block
            ym = H * xm;
            // Loop over k,m idxs
            for (const auto &[km, akm] : akm_idx) {
                const auto &[k, m] = km;
                const auto &bkm = bkm_idx[km];
                int gamma_km = std::abs(k * this->Nr - m * this->Nd);
                if (akm != -1)
                    output[gamma_km] = {ym[akm], ym[bkm]};
            }
        }
        return output;
    }

    double get_wo_dot() const { return wo_dot; }

    double get_df() const { return this->df; }

    void compute_normal_matrix() override final {
        // Compute weighted normal matrix
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            this->N.at(m0) = H.at(m0).transpose() *
                             (1 / this->n_rgt * Pyy.at(m0))
                                 .diagonal()
                                 .cwiseInverse()
                                 .asDiagonal() *
                             H.at(m0);
        }
    };
};

#endif //_BASE_OBSERVATION_HPP_