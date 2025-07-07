
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <great>

TEST(TestConstellation, BuildMatrix) {
    int Nr = 31;
    int Nd = 2;
    int Np = 4;
    int l_max = Nr / 2 * Np;
    Eigen::VectorXd inclinations(Np);
    inclinations << 90, 80, 70, 60;
    inclinations = inclinations * M_PI / 180;
    Constellation constellation(l_max, inclinations, 200e3, Nr, Nd);
    auto range_psd = [](double f) {
        return 2.62 * sqrt(1 + pow(0.003 / f, 2)) * 1e-6;
    };
    auto acc_psd = [](double f) {
        return 1e-11 * sqrt((1e-1 / f) + 1 + pow(f / 1.5e-1, 4));
    };

    constellation.set_observation_error(range_psd, acc_psd);
    constellation.set_kaula_regularization();
    constellation.block_solve();
    auto functional = EquivalentWaterHeight(0);
    double sigma = constellation.synthesis_average(functional);

    std::cout << sigma << std::endl;
}
