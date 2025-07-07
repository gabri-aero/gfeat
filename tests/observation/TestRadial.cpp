#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"

#include <gfeat>

double deg2rad(double deg) { return deg * M_PI / 180.0; }

TEST(TestRadial, BuildMatrix) {
    int Nr = 16;
    int Nd = 1;
    double I = deg2rad(30.082176279599226895);
    int l_max = 4;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Radial radial(l_max, I, Nr, Nd, deg2rad(0.06990205365958082219),
                  deg2rad(0.0014901375793001010505));
    radial.set_observation_error(allan_variance);
    GravityField gravity_field("GOCO05c.gfc", l_max, l_max);
    radial.simulate_observations(gravity_field);
    Eigen::VectorXd x(3);
    x << -4.841694588430318E-04, +9.571897944720510E-07, +5.399926430982187E-07;
    std::cout << radial.get_H_blocks()[0] << std::endl;
    std::cout << radial.get_H_blocks()[0] * x << std::endl;
}
