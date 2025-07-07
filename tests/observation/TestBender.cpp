#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include <gfeat>
#include <gtest/gtest.h>

TEST(TestBender, BuildMatrix) {
    Eigen::VectorXd I(2);
    I << M_PI_2, 1;
    double eta = 100 * M_PI / 180;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Bender bender(30, I, eta, 15, 1, LongitudePolicy::OVERLAPPING);
    bender.set_observation_error(allan_variance,
                                 [](double f) -> double { return 0; });
    bender.set_kaula_regularization();
    bender.block_solve();
}
