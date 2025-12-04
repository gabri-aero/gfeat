#include <chrono>
#include <gfeat>
#include <iostream>

#include "gtest/gtest.h"

TEST(TestSH, Load) {
    SphericalHarmonicsCovariance sh_covariance(20);
    sh_covariance.from_normal(
        "gravity/monthly/normals/ITSG-Grace_operational_n96_2024-12.snx");
    auto [lon, lat, sigma_ga] =
        sh_covariance.synthesis(40, 90, GravityAnomaly());
    std::cout << sigma_ga.col(0) << std::endl;
}