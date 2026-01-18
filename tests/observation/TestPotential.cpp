
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <functions>
#include <gfeat>

double deg2rad(double deg) { return deg * M_PI / 180.0; }

TEST(TestPotential, Simulate) {
    int Nr = 16;
    int Nd = 1;
    double I = deg2rad(40);
    int l_max = 2;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Potential potential(l_max, Nr, Nd, I, 0, 0);
    potential.set_observation_error(allan_variance);
    GravityField gravity_field("GOCO05c.gfc", 2, 2);
    potential.simulate_observations(gravity_field);

    Flmp flmp(2, I);
    double r = potential.get_radius();
    std::cout << potential.get_H_blocks()[0] << std::endl;
    ASSERT_EQ(potential.get_H_blocks()[0](0, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmk(2, 0, 2) + flmp.get_Flmk(2, 0, -2)));
    ASSERT_EQ(potential.get_H_blocks()[0](1, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmk(2, 0, 2) - flmp.get_Flmk(2, 0, -2)));
}
