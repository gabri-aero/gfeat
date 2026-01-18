
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <gfeat>

TEST(TestBender, Bender) {
    logger.set_verbosity(Verbosity::Info);
    // Define settings for Bender
    int l_max = 100;
    int Nr = 78;
    int Nd = 5;
    double rho_0 = 290e3;
    Eigen::VectorX<double> I(2);
    I << M_PI_2, M_PI / 180 * 60;
    // Define instrument noise ASD
    std::function<double(double)> lri = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    std::function<double(double)> acc = [](double f) -> double {
        return 2e-12 * sqrt((2e-3 / f) + 1 + pow(f / 1e-1, 4));
    };
    // Build observation system and solve
    Constellation bender(l_max, Nr, Nd, I, rho_0,
                         LongitudePolicy::INTERLEAVING);
    bender.set_observation_error(lri, acc);
    bender.set_kaula_regularization();
    bender.block_solve();
    // Compute EWH RMS on the sphere
    std::cout << "EWH RMS: "
              << bender.synthesis_average(EquivalentWaterHeight(0.0));
}
