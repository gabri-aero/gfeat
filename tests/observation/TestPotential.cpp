
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <functions>
#include <great>

double deg2rad(double deg) { return deg * M_PI / 180.0; }

TEST(TestPotential, BuildMatrix) {
    int Nr = 16;
    int Nd = 1;
    double I = deg2rad(40);
    int l_max = 2;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Potential potential(l_max, I, Nr, Nd, 0, 0);
    potential.populate_design_matrix();
    potential.set_observation_error(allan_variance);
    GravityField gravity_field("GOCO05c.gfc", 2, 2);
    potential.simulate_observations(gravity_field);

    Flmp flmp(2, I);
    double r = potential.get_radius();
    ASSERT_EQ(potential.get_H_blocks()[0](0, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmp(2, 0, 0) + flmp.get_Flmp(2, 0, 2)));
    ASSERT_EQ(potential.get_H_blocks()[0](1, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmp(2, 0, 0) - flmp.get_Flmp(2, 0, 2)));

    ASSERT_EQ(potential.get_H_blocks()[1](0, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmp(2, 0, 0) + flmp.get_Flmp(2, 0, 2)));
    ASSERT_EQ(potential.get_H_blocks()[1](1, 0),
              mu / r * pow(ae / r, 2) *
                  (flmp.get_Flmp(2, 0, 0) - flmp.get_Flmp(2, 0, 2)));
}
