#include <gfeat>
#include <iostream>

#include "gtest/gtest.h"

TEST(TestGravityField, LoadGravityField) {
    auto sh = GravityField(
        std::string(DATA_DIR) + "/gravity/static/ITSG-Grace2018s.gfc", 100);
    std::cout << sh.Clm(2, 0) << std::endl;
    auto sh2 = GravityField(
        std::string(DATA_DIR) + "/gravity/static/GOCO05c.gfc", 100);
    std::cout << sh2.Clm(2, 0) << std::endl;
}

TEST(TestGravityField, Potential) {
    int l_max = 90;
    GravityField gravity_field(
        std::string(DATA_DIR) + "/gravity/static/GOCO05c.gfc", l_max, l_max);
    Eigen::Vector3d r_ecrf(7000e3, 1000e3, 0);
    std::cout << gravity_field.potential(r_ecrf) << std::endl;
    std::cout << gravity_field.mu << std::endl;
    std::cout << gravity_field.ae << std::endl;
    std::cout << 3.984e14 / 7000e3 * (0.00048) << std::endl;
}