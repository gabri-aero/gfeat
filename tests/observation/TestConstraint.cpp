#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"

#include <gfeat>

double deg2rad(double deg) { return deg * M_PI / 180.0; }

TEST(TestRadial, DesignMatrix) {
    logger.set_verbosity(Verbosity::Info);
    double I = M_PI_2;
    double rho_0 = 200e3;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Collinear collinear(25, 16, 1, I, rho_0, 0.1);
    collinear.set_observation_error(allan_variance);
    collinear.set_kaula_regularization();
    collinear.set_Clm_constraint(2, 0, 0.2e-10);
    collinear.block_solve();
    std::cout << "Sigma C20: " << collinear.get_sigma_Clm(2, 0);
}
