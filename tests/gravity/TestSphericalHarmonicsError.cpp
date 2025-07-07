#include <chrono>
#include <great>
#include <iostream>

#include "gtest/gtest.h"

TEST(TestSH, Load) {
    SphericalHarmonicsError sh_error(
        20, "gravity/monthly/normals/ITSG-Grace_operational_n96_2024-12.snx");
    auto [lon, lat, sigma_ga] = sh_error.synthesis(40, 90, GravityAnomaly());
    std::cout << sigma_ga.col(0) << std::endl;
}