#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include "../Planet.hpp"

#include <Eigen/Dense>
#include <iostream>

Eigen::Matrix<double, 3, 3> R1(double a) {
    return Eigen::Matrix<double, 3, 3>{{1, 0, 0},
                                       {0, std::cos(a), std::sin(a)},
                                       {0, -std::sin(a), std::cos(a)}};
}

Eigen::Matrix<double, 3, 3> R2(double a) {
    return Eigen::Matrix<double, 3, 3>{{std::cos(a), 0, -std::sin(a)},
                                       {0, 1, 0},
                                       {std::sin(a), 0, std::cos(a)}};
}

Eigen::Matrix<double, 3, 3> R3(double a) {
    return Eigen::Matrix<double, 3, 3>{{std::cos(a), std::sin(a), 0},
                                       {-std::sin(a), std::cos(a), 0},
                                       {0, 0, 1}};
}

Eigen::Vector3d cart2sph(Eigen::Vector3d r_cart) {
    // Retrieve vector components
    double x = r_cart(0);
    double y = r_cart(1);
    double z = r_cart(2);
    // Compute norm
    double r = r_cart.norm();
    // Compute longitude
    double lon = atan2(y, x);
    // Compute latitude
    double lat = atan2(z, sqrt(x * x + y * y));
    // Return
    return Eigen::Vector3d{r, lat, lon};
}

Eigen::Vector3d sph2cart(Eigen::Vector3d sph) {
    // Retrieve vector components
    double r = sph(0);
    double lon = sph(2);
    double lat = sph(1);
    // Define cartesian position
    Eigen::Vector3d r_cart;
    r_cart(0) = r * cos(lat) * cos(lon);
    r_cart(1) = r * cos(lat) * sin(lon);
    r_cart(2) = r * sin(lat);
    // Return
    return r_cart;
}

Eigen::Matrix3d dsph_dcart(Eigen::Vector3d cart) {
    // Pre-allocate
    Eigen::Matrix3d dsph_dcart;
    // Unary vectors
    Eigen::Vector3d ux{1, 0, 0};
    Eigen::Vector3d uy{0, 1, 0};
    Eigen::Vector3d uz{0, 0, 1};
    // Compute derivative vectors
    double x2_y2 = pow(cart(0), 2) + pow(cart(1), 2);
    double r = cart.norm();
    dsph_dcart.col(0) = cart / r; // dr_dcart
    dsph_dcart.col(1) =
        (uz - cart * cart(2) / pow(r, 2)) / sqrt(x2_y2);       // dphi_dcart
    dsph_dcart.col(2) = (cart(0) * uy - cart(1) * ux) / x2_y2; // dlam_dcart
    return dsph_dcart;
}

Eigen::Matrix3d dcart_dsph(Eigen::Vector3d cart) {
    return dsph_dcart(cart).transpose();
}

std::array<Eigen::Matrix3d, 3> ddsph_ddcart(Eigen::Vector3d cart) {
    std::array<Eigen::Matrix3d, 3> ddsph_ddcart;
    // Define constants
    double x = cart(0);
    double y = cart(1);
    double z = cart(2);
    double r = cart.norm();
    double r2 = pow(r, 2);
    double x2 = pow(x, 2);
    double y2 = pow(y, 2);
    double z2 = pow(z, 2);
    double r3 = pow(r, 3);
    double r4 = pow(r, 4);
    double rho = sqrt(x2 + y2);
    double rho3 = pow(rho, 3);
    double rho2 = pow(rho, 2);
    double rho4 = pow(rho, 4);
    // Radial components
    Eigen::Matrix3d &ddr_ddcart = ddsph_ddcart[0];
    ddr_ddcart(0, 0) = (r2 - x2) / r3;
    ddr_ddcart(0, 1) = -x * y / r3;
    ddr_ddcart(0, 2) = -x * z / r3;
    ddr_ddcart(1, 0) = ddr_ddcart(0, 1);
    ddr_ddcart(1, 1) = (r2 - y2) / r3;
    ddr_ddcart(1, 2) = -y * z / r3;
    ddr_ddcart(2, 0) = ddr_ddcart(0, 2);
    ddr_ddcart(2, 1) = ddr_ddcart(1, 2);
    ddr_ddcart(2, 2) = (r2 - z2) / r3;
    // Latitude components
    Eigen::Matrix3d &ddphi_ddcart = ddsph_ddcart[1];
    ddphi_ddcart(0, 0) =
        z * (2 * x2 * x2 + x2 * y2 - y2 * (y2 + z2)) / (rho3 * r4);
    ddphi_ddcart(0, 1) = x * y * z * (3 * rho2 + z2) / (rho3 * r4);
    ddphi_ddcart(0, 2) = -x * (rho2 - z2) / (rho * r4);
    ddphi_ddcart(1, 0) = ddphi_ddcart(0, 1);
    ddphi_ddcart(1, 1) =
        -z * (x2 * x2 - 2 * y2 * y2 + x2 * (-y2 + z2)) / (rho3 * r4);
    ddphi_ddcart(1, 2) = -y * (rho2 - z2) / (rho * r4);
    ddphi_ddcart(2, 0) = ddphi_ddcart(0, 2);
    ddphi_ddcart(2, 1) = ddphi_ddcart(1, 2);
    ddphi_ddcart(2, 2) = -2 * rho * z / r4;
    // Longitude components
    Eigen::Matrix3d &ddlam_ddcart = ddsph_ddcart[2];
    ddlam_ddcart.setZero();
    ddlam_ddcart(0, 0) = 2 * x * y / rho4;
    ddlam_ddcart(0, 1) = (-x2 + y2) / rho4;
    ddlam_ddcart(1, 0) = ddlam_ddcart(0, 1);
    ddlam_ddcart(1, 1) = -ddlam_ddcart(0, 0);
    return ddsph_ddcart;
}

Eigen::Matrix<double, 3, 3> Omega() {
    Eigen::Matrix<double, 3, 3> S;
    S << 0, -theta_dot, 0, theta_dot, 0, 0, 0, 0, 0;
    return S;
}

Eigen::Vector3d eci2ecrf(Eigen::Vector3d r_eci, double t) {
    auto M_RI = R3(theta_dot * t);
    return M_RI * r_eci;
}

std::pair<Eigen::Vector3d, Eigen::Vector3d>
eci2ecrf(Eigen::Vector3d r_eci, Eigen::Vector3d v_eci, double t) {
    auto M_RI = R3(theta_dot * t);
    Eigen::Matrix<double, 3, 3> S_omega = Omega();
    auto r_ecrf = eci2ecrf(r_eci, t);
    auto v_ecrf = M_RI * v_eci - S_omega * r_ecrf;
    return {r_ecrf, v_ecrf};
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
eci2ecrf(Eigen::Vector3d r_eci, Eigen::Vector3d v_eci, Eigen::Vector3d a_eci,
         double t) {
    auto M_RI = R3(theta_dot * t);
    Eigen::Matrix<double, 3, 3> S_omega = Omega();
    auto [r_ecrf, v_ecrf] = eci2ecrf(r_eci, v_eci, t);
    auto a_ecrf =
        M_RI * a_eci - 2 * S_omega * v_ecrf - S_omega * S_omega * r_ecrf;
    return {r_ecrf, v_ecrf, a_ecrf};
}

Eigen::Vector3d ecrf2eci(Eigen::Vector3d r_ecrf, double t) {
    auto M_IR = R3(-theta_dot * t);
    return M_IR * r_ecrf;
}

std::pair<Eigen::Vector3d, Eigen::Vector3d>
ecrf2eci(Eigen::Vector3d r_ecrf, Eigen::Vector3d v_ecrf, double t) {
    auto M_IR = R3(-theta_dot * t);
    Eigen::Matrix<double, 3, 3> S_omega = Omega();
    auto r_eci = ecrf2eci(r_ecrf, t);
    auto v_eci = M_IR * (v_ecrf + S_omega * r_ecrf);
    return {r_eci, v_eci};
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
ecrf2eci(Eigen::Vector3d r_ecrf, Eigen::Vector3d v_ecrf, Eigen::Vector3d a_ecrf,
         double t) {
    auto M_IR = R3(-theta_dot * t);
    Eigen::Matrix<double, 3, 3> S_omega = Omega();
    auto [r_eci, v_eci] = ecrf2eci(r_ecrf, v_ecrf, t);
    auto a_eci =
        M_IR * (a_ecrf + 2 * S_omega * v_ecrf + S_omega * S_omega * r_ecrf);
    return {r_eci, v_eci, a_eci};
}

Eigen::Matrix<double, 6, 6> eci2ecrf(Eigen::Matrix<double, 6, 6> phi_eci,
                                     double t) {
    Eigen::Matrix<double, 6, 6> dx_ecrf_dx_eci_t, dx_eci_dx_ecrf_0;
    // Initialize to zero
    dx_ecrf_dx_eci_t.setZero();
    dx_eci_dx_ecrf_0.setZero();
    // Compute necessary sub-matrices
    auto I = Eigen::Matrix<double, 3, 3>::Identity();
    auto M_RI = R3(theta_dot * t);
    Eigen::Matrix<double, 3, 3> S_omega = Omega();
    // Assign submatrices
    dx_eci_dx_ecrf_0.block(0, 0, 3, 3) = I;
    dx_eci_dx_ecrf_0.block(3, 3, 3, 3) = I;
    dx_eci_dx_ecrf_0.block(3, 0, 3, 3) = S_omega;
    // Assign submatrices
    dx_ecrf_dx_eci_t.block(0, 0, 3, 3) = M_RI;
    dx_ecrf_dx_eci_t.block(3, 3, 3, 3) = M_RI;
    dx_ecrf_dx_eci_t.block(3, 0, 3, 3) = -S_omega * M_RI;
    // Apply transformation to STM
    return dx_ecrf_dx_eci_t * phi_eci * dx_eci_dx_ecrf_0;
}

#endif // _COORDINATES_HPP_