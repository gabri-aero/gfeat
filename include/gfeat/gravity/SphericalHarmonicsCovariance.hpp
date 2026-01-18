#ifndef _SPHERICAL_HARMONICS_COVARIANCE_HPP_
#define _SPHERICAL_HARMONICS_COVARIANCE_HPP_

#define EIGEN_HAS_OPENMP // Force Eigen to see OpenMP
#include <omp.h>

#include "../utils/Logger.hpp"
#include "../utils/Sinex.hpp"
#include "Functionals.hpp"

#include <chrono>
#include <fft>
#include <functions>
#include <tuple>

class SphericalHarmonicsCovariance {
    int l_max;

    int clm_idx(int l, int m) const {
        return m == 0 ? l * l - 4 : l * l + 2 * m - 1 - 4;
    }

    int slm_idx(int l, int m) const { return l * l + 2 * m - 4; }

public:
    Eigen::MatrixXd Pxx;
    SphericalHarmonicsCovariance(int l_max) : l_max(l_max) {
        int n = clm_idx(l_max + 1, 0);
        this->Pxx.resize(n, n);
    }
    SphericalHarmonicsCovariance from_normal(std::string filename,
                                             double scaling_factor = 1.0) {
        auto start = std::chrono::high_resolution_clock::now();
        Sinex sinex(filename);
        auto N = sinex.get_normal_matrix(l_max);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> s = end - start;
        logger << "Normal matrix loaded in " << s.count() << " secs \n";
        start = std::chrono::high_resolution_clock::now();
        Eigen::LDLT<Eigen::MatrixXd> ldlt(N);
        this->Pxx = Eigen::MatrixXd::Identity(N.rows(), N.cols());
        ldlt.solveInPlace(this->Pxx);
        this->Pxx *= scaling_factor;
        end = std::chrono::high_resolution_clock::now();
        s = end - start;
        logger << "Normal matrix inverted in " << s.count() << " secs \n";
        return *this;
    }

    Eigen::MatrixXd get_Pxx() const { return Pxx; }

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
    synthesis(int n_lon, int n_lat,
              const BaseFunctional &gravity_field_functional) {
        auto start = std::chrono::high_resolution_clock::now();
        // Adjust number of longitude nodes for FFT using full SH spectrum
        if (n_lon / 2 < l_max) {
            n_lon = 2 * l_max;
            std::cerr << "Warning: Number of longitude nodes increased to "
                      << n_lon
                      << " reach full SH spectrum during FFT synthesis over "
                         "parallels."
                      << std::endl;
        }
        // Create longitude and latitude vectors
        double dlon = 2 * M_PI / n_lon;
        double dlat = M_PI / n_lat;
        Eigen::VectorXd lon_vec =
            Eigen::VectorXd::LinSpaced(n_lon, 0, 2 * M_PI - dlon);
        auto lon_vec_wrapped =
            Eigen::VectorXd::LinSpaced(n_lon, -M_PI, M_PI - dlon);
        auto lat_vec = Eigen::VectorXd::LinSpaced(n_lat, -M_PI / 2 + dlat / 2,
                                                  M_PI / 2 - dlat / 2);
        // Apply broadcasting for generating mesh
        Eigen::MatrixXd lon_mesh =
            (lon_vec_wrapped * 180 / M_PI).transpose().replicate(n_lat, 1);
        Eigen::MatrixXd lat_mesh = (lat_vec * 180 / M_PI).replicate(1, n_lon);
        // Pre-allocate vector of partials and W matrix
        int nx = Pxx.rows();
        Eigen::MatrixXd W(nx, n_lon);
        Eigen::VectorXd v(nx);
        v.setZero();
        W.setZero();
        // Define degree term
        Eigen::VectorXd f(l_max + 1);
        for (int l = 0; l <= l_max; l++) {
            f(l) = gravity_field_functional.degree_dependent_factor(l);
        }
        f = f * gravity_field_functional.common_factor;
        // Pre-allocate cos(m lam), sin(m lam) terms
        Eigen::MatrixXd cos_mlam(l_max + 1, n_lon);
        Eigen::MatrixXd sin_mlam(l_max + 1, n_lon);
        for (int m = 0; m <= l_max; m++) {
            for (int lon_idx = 0; lon_idx < n_lon; lon_idx++) {
                cos_mlam(m, lon_idx) = cos(m * lon_vec(lon_idx));
                sin_mlam(m, lon_idx) = sin(m * lon_vec(lon_idx));
            }
        }
        // Null out C20 components
        Pxx.row(0).setZero();
        Pxx.col(0).setZero();
        // Pre-allocate functional std deviation
        Eigen::MatrixXd y(n_lat, n_lon);
        // Define Am, Bm
        int n = n_lon / 2;
        Eigen::VectorXd A(n + 1), B(n + 1);
        // Define latitude loop
        logger << "Covariance propagation progress" << std::endl;
        for (int lat_idx = 0; lat_idx < n_lat; lat_idx++) {
            logger << lat_idx << "/" << n_lat << "\r" << std::flush;
            // Compute Legendre polynomials
            Plm plm(l_max, M_PI_2 - lat_vec(lat_idx));
            // Loop over covariance rows
            for (int row_idx = 0; row_idx < nx; row_idx++) {
                // Fill A, B (latitude dependent terms)
                A.setZero();
                B.setZero();
                for (int m = 0; m <= l_max; m++) {
                    for (int l = std::max(m, 2); l <= l_max; l++) {
                        A(m) += f(l) * plm.get_Plm_bar(l, m) *
                                Pxx(row_idx, clm_idx(l, m));
                        if (m > 0)
                            B(m) += f(l) * plm.get_Plm_bar(l, m) *
                                    Pxx(row_idx, slm_idx(l, m));
                    }
                }
                // Build FFT complex vector
                Eigen::VectorX<std::complex<double>> fhat(n_lon);
                fhat.setZero();
                fhat(0) = A(0);
                auto j = std::complex<double>(0.0, 1.0);
                for (int i = 1; i <= n; i++) {
                    fhat(i) = fhat(i) + 0.5 * (A(i) - j * B(i));
                    fhat(n_lon - i) = fhat(n_lon - i) + 0.5 * (A(i) + j * B(i));
                }
                // Inverse FFT call to compute all Wi longitude terms at once
                fhat *= n_lon;
                auto wi = ifft(fhat);
                std::transform(
                    wi.begin(), wi.end(), W.row(row_idx).begin(),
                    [](const std::complex<double> &z) { return z.real(); });
            }
            // Loop over longitudes
            for (int lon_idx = 0; lon_idx < n_lon; lon_idx++) {
                // Fill vector of partials
                for (int l = 2; l <= l_max; l++) {
                    for (int m = 0; m <= l_max; m++) {
                        v(clm_idx(l, m)) =
                            f(l) * plm.get_Plm_bar(l, m) * cos_mlam(m, lon_idx);
                        if (m > 0)
                            v(slm_idx(l, m)) = f(l) * plm.get_Plm_bar(l, m) *
                                               sin_mlam(m, lon_idx);
                    }
                }
                // Compute std deviation
                y(lat_idx, lon_idx) = sqrt(v.transpose() * W.col(lon_idx));
            }
            // Rotate to match the [-pi, pi] range
            std::rotate(y.row(lat_idx).begin(),
                        y.row(lat_idx).begin() + n_lon / 2,
                        y.row(lat_idx).end());
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> s = end - start;
        logger << "Covariance propagation completed in " << s.count()
               << " secs \n";
        // Return values
        return {lon_mesh, lat_mesh, y};
    }

    double sigma_Clm(int l, int m) const {
        int index = clm_idx(l, m);
        return this->Pxx(index, index);
    }

    double sigma_Slm(int l, int m) const {
        int index = slm_idx(l, m);
        return this->Pxx(index, index);
    }

    Eigen::VectorXd degree_variance() const {
        Eigen::VectorXd degree_variance(this->l_max + 1);
        degree_variance.setZero();
        for (int l = 2; l <= l_max; l++) {
            for (int m = 0; m <= l; m++) {
                if (m == 0) {
                    degree_variance(l) += sigma_Clm(l, m);
                } else {
                    degree_variance(l) += sigma_Clm(l, m) + sigma_Slm(l, m);
                }
            }
        }
        return degree_variance;
    }

    Eigen::VectorXd rms_per_coefficient_per_degree() const {
        Eigen::VectorXd result(l_max + 1);
        auto deg_var = this->degree_variance();
        for (int l = 2; l <= this->l_max; l++) {
            result(l) = std::sqrt(deg_var.coeff(l) / (2 * l + 1));
        }
        return result;
    }
};

#endif // SPHERICAL_HARMONICS_COVARIANCE_HPP_