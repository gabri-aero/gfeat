#ifndef _SPHERICAL_HARMONICS_HPP_
#define _SPHERICAL_HARMONICS_HPP_

#include "../utils/Coordinates.hpp"
#include "../utils/Logger.hpp"
#include "Functionals.hpp"

#include <Eigen/Dense>
#include <fft>
#include <functional>
#include <functions>
#include <iostream>

class SphericalHarmonics {
  public:
    // Kronecker delta definition
    int d(int i, int j) { return i == j ? 1 : 0; }

    // Auxiliar function for general purpose scalar observations
    double aux_function(Eigen::Vector3d r_ecrf,
                        std::function<double(int, double)> f,
                        std::function<double(double)> g) {
        // Convert r_ecrf to polar coordinates
        auto polar = cart2sph(r_ecrf);
        // Retrieve polar coomponents
        double r = polar(0);
        double lam = polar(2);       // longitude
        double phi = polar(1);       // latitude
        double theta = M_PI_2 - phi; // colatitude
        // Compute associated Legendre polynomials
        Plm plm{l_max, theta, false};
        // Pre-compute some values
        Eigen::VectorXd cos_mlam(l_max + 1);
        Eigen::VectorXd sin_mlam(l_max + 1);
        for (int m = 0; m <= l_max; m++) {
            cos_mlam(m) = cos(m * lam);
            sin_mlam(m) = sin(m * lam);
        }
        // Allocate variables
        double V = 0;    // Potential
        double Vl = 0;   // Potential common term over degree
        double clm, slm; // Stokes coefficients
        // Loop over degree
        for (int l = 2; l <= l_max; l++) {
            // Reset degree common term
            Vl = 0;
            for (int m = 0; m <= std::min(l, m_max); m++) {
                // Get Stokes coefficients
                clm = Clm(l, m);
                slm = m > 0 ? Slm(l, m) : 0;
                // Compute order contribution to degree common term
                Vl += plm.get_Plm_bar(l, m) *
                      (clm * cos_mlam(m) + slm * sin_mlam(m));
            }
            // Add degree constant
            Vl *= f(l, r);
            // Append contribution to total potential
            V += Vl;
        }
        // Add global factor
        V *= g(r);
        return V;
    };

    Eigen::Vector3d dV_dsph(Eigen::Vector3d cart) {
        // Convert r_ecrf to polar coordinates
        auto polar = cart2sph(cart);
        // Retrieve polar coomponents
        double r = polar(0);
        double lam = polar(2);       // longitude
        double phi = polar(1);       // latitude
        double theta = M_PI_2 - phi; // colatitude
        // Pre-allocate variables
        double dU_dr = 0;
        double dU_dphi = 0;
        double dU_dlam = 0;
        double dU_dr_l, dU_dphi_l, dU_dlam_l;
        // Compute associated Legendre polynomials
        Plm plm(l_max + 1, theta, true);
        // Pre-compute some values
        Eigen::VectorXd cos_mlam(l_max + 1);
        Eigen::VectorXd sin_mlam(l_max + 1);
        for (int m = 0; m <= l_max; m++) {
            cos_mlam(m) = cos(m * lam);
            sin_mlam(m) = sin(m * lam);
        }
        double r_hat = ae / r;
        double f;        // attenuation factor
        double clm, slm; // Stokes coefficients
        int l, m;

        // Double summation
        for (l = 2; l <= l_max; l++) { // gravity from perturbing potential
            f = pow(r_hat, l);
            // Reset variables
            dU_dr_l = 0;
            dU_dlam_l = 0;
            dU_dphi_l = 0;
            // Derivative formula (without common terms)
            for (m = 0; m <= std::min(l, m_max); m++) {
                // Get Stokes coefficients
                clm = Clm(l, m);
                slm = m > 0 ? Slm(l, m) : 0;
                // Compute contribution to partials
                dU_dr_l += plm.get_Plm_bar(l, m) *
                           (clm * cos_mlam(m) + slm * sin_mlam(m));
                dU_dphi_l -= plm.get_dPlm_bar(l, m) *
                             (clm * cos_mlam(m) + slm * sin_mlam(m));
                dU_dlam_l += m * plm.get_Plm_bar(l, m) *
                             (slm * cos_mlam(m) - clm * sin_mlam(m));
            }

            dU_dr -= (l + 1) * f * dU_dr_l / r;
            dU_dphi += f * dU_dphi_l;
            dU_dlam += f * dU_dlam_l;
        }
        dU_dr *= mu / r;
        dU_dphi *= mu / r;
        dU_dlam *= mu / r;
        return Eigen::Vector3d{dU_dr, dU_dphi, dU_dlam};
    }

    Eigen::Matrix3d ddV_ddsph(Eigen::Vector3d cart) {
        // Convert r_ecrf to polar coordinates
        auto polar = cart2sph(cart);
        // Retrieve polar coomponents
        double r = polar(0);
        double lam = polar(2);       // longitude
        double phi = polar(1);       // latitude
        double theta = M_PI_2 - phi; // colatitude
        // Pre-allocate variables
        double ddU_ddr = 0;
        double ddU_ddphi = 0;
        double ddU_ddlam = 0;
        double ddU_dr_dphi = 0;
        double ddU_dr_dlam = 0;
        double ddU_dphi_dlam = 0;
        double ddU_ddr_l, ddU_ddphi_l, ddU_ddlam_l, ddU_dr_dphi_l,
            ddU_dr_dlam_l, ddU_dphi_dlam_l;
        // Compute associated Legendre polynomials
        Plm plm(l_max + 1, theta, true, true);
        // Pre-compute some values
        Eigen::VectorXd cos_mlam(l_max + 1);
        Eigen::VectorXd sin_mlam(l_max + 1);
        for (int m = 0; m <= l_max; m++) {
            cos_mlam(m) = cos(m * lam);
            sin_mlam(m) = sin(m * lam);
        }
        double r_hat = ae / r;
        double f;        // attenuation factor
        double clm, slm; // Stokes coefficients
        int l, m;

        // Double summation
        for (l = 2; l <= l_max; l++) { // gravity from perturbing potential
            f = pow(r_hat, l);
            // Reset variables
            ddU_ddr_l = 0;
            ddU_ddlam_l = 0;
            ddU_ddphi_l = 0;
            ddU_dr_dphi_l = 0;
            ddU_dr_dlam_l = 0;
            ddU_dphi_dlam_l = 0;
            // Derivative formula (without common terms)
            for (m = 0; m <= std::min(l, m_max); m++) {
                // Get Stokes coefficients
                clm = Clm(l, m);
                slm = m > 0 ? Slm(l, m) : 0;
                // Compute contribution to partials
                ddU_ddr_l += plm.get_Plm_bar(l, m) *
                             (clm * cos_mlam(m) + slm * sin_mlam(m));
                ddU_ddphi_l += plm.get_ddPlm_bar(l, m) *
                               (clm * cos_mlam(m) + slm * sin_mlam(m));
                ddU_ddlam_l += pow(m, 2) * plm.get_Plm_bar(l, m) *
                               (clm * cos_mlam(m) + slm * sin_mlam(m));
                ddU_dr_dphi_l += plm.get_dPlm_bar(l, m) *
                                 (clm * cos_mlam(m) + slm * sin_mlam(m));
                ddU_dr_dlam_l += plm.get_Plm_bar(l, m) * m *
                                 (slm * cos_mlam(m) - clm * sin_mlam(m));
                ddU_dphi_dlam_l += plm.get_dPlm_bar(l, m) * m *
                                   (slm * cos_mlam(m) - clm * sin_mlam(m));
            }
            ddU_ddr += (l + 1) * (l + 2) * f * ddU_ddr_l;
            ddU_ddphi += f * ddU_ddphi_l;
            ddU_ddlam -= f * ddU_ddlam_l;
            ddU_dr_dphi += (l + 1) * f * ddU_dr_dphi_l;
            ddU_dr_dlam -= (l + 1) * f * ddU_dr_dlam_l;
            ddU_dphi_dlam -= f * ddU_dphi_dlam_l;
        }
        // Add common factor
        ddU_ddr *= mu / pow(r, 3);
        ddU_ddphi *= mu / r;
        ddU_ddlam *= mu / r;
        ddU_dr_dphi *= mu / pow(r, 2);
        ddU_dr_dlam *= mu / pow(r, 2);
        ddU_dphi_dlam *= mu / r;
        // Build full tensor
        Eigen::Matrix3d mat;
        mat << ddU_ddr, ddU_dr_dphi, ddU_dr_dlam, ddU_dr_dphi, ddU_ddphi,
            ddU_dphi_dlam, ddU_dr_dlam, ddU_dphi_dlam, ddU_ddlam;
        return mat;
    }

  public:
    // Class public attributes
    int l_max;
    int m_max;
    double mu;
    double ae;
    Eigen::MatrixXd coefficients;

    // Constructors
    SphericalHarmonics() = default;
    SphericalHarmonics(int l_max) : SphericalHarmonics(l_max, l_max) {}
    SphericalHarmonics(int l_max, int m_max) {
        this->l_max = l_max;
        this->m_max = m_max;
        // Allocate coefficient matrices
        this->coefficients.resize(l_max + 1, l_max + 1);
        // Assign default
        this->coefficients.setZero();
        this->Clm(0, 0) = 1; // by definition
    }

    // Stokes coefficients indexing functions
    double Clm(int l, int m) const { return coefficients(l, m); }
    double Slm(int l, int m) const { return coefficients(m - 1, l); }
    double &Clm(int l, int m) { return coefficients(l, m); }
    double &Slm(int l, int m) { return coefficients(m - 1, l); }

    Eigen::Vector3d gravity(Eigen::Vector3d r_ecrf) {
        return dcart_dsph(r_ecrf).transpose() *
               dV_dsph(r_ecrf); // Gradient covariant transform
    }

    Eigen::Matrix3d gravity_gradient(Eigen::Vector3d r_ecrf) {
        Eigen::Matrix3d partial_term = dsph_dcart(r_ecrf) * ddV_ddsph(r_ecrf) *
                                       dsph_dcart(r_ecrf).transpose();
        auto connect = [](const auto &a, const auto &b) {
            return a[0] * b(0) + a[1] * b(1) + a[2] * b(2);
        };
        Eigen::Matrix3d connection_term =
            connect(ddsph_ddcart(r_ecrf), dV_dsph(r_ecrf));
        return partial_term + connection_term;
    }

    double potential(Eigen::Vector3d r_ecrf) {
        return aux_function(
            r_ecrf,
            [this](int l, double r) -> double { return pow(ae / r, l); },
            [this](double r) -> double { return mu / r; });
    }

    double geoid_height(Eigen::Vector3d r_ecrf) {
        return aux_function(
            r_ecrf, [](int l, double r) -> double { return 1; },
            [this](double r) -> double { return ae; });
    }

    double gravity_anomaly(Eigen::Vector3d r_ecrf) {
        return aux_function(
            r_ecrf, [](int l, double r) -> double { return (l - 1); },
            [this](double r) -> double { return 1e5 * mu / pow(ae, 2); });
    }

    Eigen::VectorXd degree_variance() const {
        Eigen::VectorXd degree_variance(this->l_max + 1);
        degree_variance.setZero();
        for (int l = 2; l <= l_max; l++) {
            for (int m = 0; m <= l; m++) {
                if (m == 0) {
                    degree_variance(l) += pow(Clm(l, m), 2);
                } else {
                    degree_variance(l) += pow(Clm(l, m), 2) + pow(Slm(l, m), 2);
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

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
    synthesis(int n_lon, int n_lat,
              const BaseFunctional &gravity_field_functional) {
        // Adjust number of longitude nodes for FFT using full SH spectrum
        if (n_lon / 2 < l_max) {
            n_lon = 2 * l_max;
            std::cerr << "Warning: Number of longitude nodes increased to "
                      << 2 * l_max
                      << " reach full SH spectrum during FFT synthesis over "
                         "parallels."
                      << std::endl;
        }
        // Create longitude and latitude vectors
        auto lon_vec = Eigen::VectorXd::LinSpaced(n_lon, -M_PI, M_PI);
        auto lat_vec = Eigen::VectorXd::LinSpaced(n_lat, -M_PI / 2, M_PI / 2);
        // Apply broadcasting for generating mesh
        Eigen::MatrixXd lon_mesh =
            (lon_vec * 180 / M_PI).transpose().replicate(n_lat, 1);
        Eigen::MatrixXd lat_mesh = (lat_vec * 180 / M_PI).replicate(1, n_lon);
        // Define observation
        Eigen::MatrixXd y(n_lat, n_lon);
        y.setZero();
        // Define degree term
        Eigen::VectorXd f(l_max + 1);
        for (int l = 0; l <= l_max; l++) {
            f(l) = gravity_field_functional.degree_dependent_factor(l);
        }
        f = f * gravity_field_functional.common_factor;
        // Force C20 = 0
        double C20 = Clm(2, 0);
        Clm(2, 0) = 0;
        // Define Am, Bm
        int n = n_lon / 2;
        Eigen::VectorXd A(n + 1), B(n + 1);
        // Define latitude loop
        for (int lat_idx = 0; lat_idx < n_lat; lat_idx++) {
            Logger::instance() << "\rSynthesis progress: " << lat_idx + 1 << "/"
                               << n_lat << std::flush;
            // Compute Legendre polynomials
            Plm plm(l_max, M_PI_2 - lat_vec(lat_idx));
            // Fill A, B
            A.setZero();
            B.setZero();
            for (int m = 0; m <= l_max; m++) {
                for (int l = std::max(m, 2); l <= l_max; l++) {
                    A(m) += f(l) * plm.get_Plm_bar(l, m) * Clm(l, m);
                    if (m > 0)
                        B(m) += f(l) * plm.get_Plm_bar(l, m) * Slm(l, m);
                }
            }
            // Build FFT complex vector
            Eigen::VectorX<std::complex<double>> f(n_lon);
            f.setZero();
            f(0) = A(0);
            auto j = std::complex<double>(0.0, 1.0);
            for (int i = 1; i <= n; i++) {
                f(i) = f(i) + 0.5 * (A(i) - j * B(i));
                f(n_lon - i) = f(n_lon - i) + 0.5 * (A(i) + j * B(i));
            }
            // Inverse FFT call to synthesize all terms along a same parallel
            f *= n_lon;
            auto y_hat = ifft(f);
            std::transform(
                y_hat.begin(), y_hat.end(), y.row(lat_idx).begin(),
                [](const std::complex<double> &z) { return z.real(); });
            std::rotate(y.row(lat_idx).begin(),
                        y.row(lat_idx).begin() + n_lon / 2,
                        y.row(lat_idx).end());
        }
        Clm(2, 0) = C20;
        // Return observation
        return {lon_mesh, lat_mesh, y};
    }
};

#endif // _SPHERICAL_HARMONICS_HPP_
