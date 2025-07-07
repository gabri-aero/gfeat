#include <gfeat>

#include "gtest/gtest.h"

TEST(TestSH, GravityTest) {
    auto sh = GravityField("gravity/static/GOCO05c.gfc", 100);
    Eigen::Vector3d sph = {7000e3, 1, -1};
    Eigen::Vector3d r_ecrf = sph2cart(sph);
    Eigen::Vector3d r_ecrf_shifted;
    Eigen::Vector3d g = sh.gravity(r_ecrf);

    double h = 1e-4;
    for (int i = 0; i < 3; i++) {
        r_ecrf_shifted = r_ecrf;
        r_ecrf_shifted(i) += h;
        double gi = (sh.potential(r_ecrf_shifted) - sh.potential(r_ecrf)) / h;
        ASSERT_NEAR((g(i) - gi) / gi, 0, 1e-4);
    }
}

TEST(TestSH, GravityGradientSphericalTest) {
    auto sh = GravityField("gravity/static/GOCO05c.gfc", 100);
    Eigen::Vector3d sph = {7000e3, 1, -1};
    Eigen::Vector3d r_ecrf = sph2cart(sph);
    Eigen::Vector3d r_ecrf_shifted, sph_shifted;
    Eigen::Matrix3d gg = sh.ddV_ddsph(r_ecrf);
    Eigen::Matrix3d gg_num;

    Eigen::Vector3d h = {1e-5, 1e-10, 1e-10};
    for (int i = 0; i < 3; i++) {
        sph_shifted = sph;
        sph_shifted(i) += h(i);
        r_ecrf_shifted = sph2cart(sph_shifted);
        gg_num.row(i) =
            (sh.dV_dsph(r_ecrf_shifted) - sh.dV_dsph(r_ecrf)) / h(i);
        for (int j = 0; j < 3; j++) {
            ASSERT_NEAR((gg_num(i, j) - gg(i, j)) / gg(i, j), 0, 5e-4);
        }
    }
}

TEST(TestSH, GravityGradientCartesianTest) {
    auto sh = GravityField("gravity/static/GOCO05c.gfc", 100);
    Eigen::Vector3d sph = {7000e3, 1, -1};
    Eigen::Vector3d r_ecrf = sph2cart(sph);
    Eigen::Vector3d r_ecrf_shifted, sph_shifted;
    Eigen::Matrix3d gg = sh.gravity_gradient(r_ecrf);
    Eigen::Matrix3d gg_num;

    std::cout << gg << std::endl << std::endl;

    double h = 1e-3;
    for (int i = 0; i < 3; i++) {
        r_ecrf_shifted = r_ecrf;
        r_ecrf_shifted(i) += h;
        gg_num.row(i) = (sh.gravity(r_ecrf_shifted) - sh.gravity(r_ecrf)) / h;
        for (int j = 0; j < 3; j++) {
            ASSERT_NEAR((gg_num(i, j) - gg(i, j)) / gg(i, j), 0, 1e-3);
        }
    }

    std::cout << gg_num << std::endl;
}

TEST(TestSH, Grid) {
    auto sh = GravityField("gravity/static/GOCO05c.gfc", 3);
    auto func = GravityAnomaly();
    auto [lon, lat, ga] = sh.synthesis(6, 4, func);
    std::cout << ga << std::endl;
}