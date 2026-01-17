
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <gfeat>

TEST(TestCollinear, SynthesisAverage) {
    logger.set_verbosity(Verbosity::Info);
    double I = M_PI_2;
    double rho_0 = 200e3;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Collinear collinear(25, 16, 1, I, rho_0, 0.1);
    collinear.set_observation_error(allan_variance);
    collinear.set_kaula_regularization();
    collinear.block_solve();
    auto functional = EquivalentWaterHeight(200e3);
    int n_lon = 180;
    int n_lat = 90;
    auto [lon, lat, y] = collinear.synthesis(n_lon, n_lat, functional);
    double mean_y = 0;
    double dlon = 2 * M_PI / n_lon;
    double dlat = M_PI / n_lat;
    for (int i = 0; i < n_lat; i++) {
        for (int j = 0; j < n_lon; j++) {
            mean_y +=
                pow(y(i, j), 2) * cos(lat(i, j) * M_PI / 180) * dlon * dlat;
        }
    }
    mean_y /= (4 * M_PI);
    mean_y = sqrt(mean_y);
    ASSERT_NEAR((mean_y - collinear.synthesis_average(functional)) / mean_y, 0,
                1e-5);
}
