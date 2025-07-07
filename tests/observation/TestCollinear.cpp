
#include <iomanip>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "gtest/gtest.h"
#include <great>

TEST(TestCollinear, BuildMatrix) {
    double I = M_PI_2;
    double eta = 10 * M_PI / 180;
    std::function<double(double)> allan_variance = [](double f) -> double {
        return 10e-9 * sqrt(1 + pow(10e-3 / f, 2));
    };
    Collinear collinear(60, I, eta, 16, 1, 0.1);
    collinear.set_observation_error(allan_variance);
    collinear.set_kaula_regularization();
    collinear.block_solve();
    auto functional = GeoidHeight();
    int n_lon = 720;
    int n_lat = 360;
    auto [lon, colat, y] = collinear.synthesis(n_lon, n_lat, functional);
    double mean_y = 0;
    double dlon = 2 * M_PI / n_lon;
    double dlat = M_PI / n_lat;
    for (int i = 0; i < n_lat; i++) {
        Plm plm(2, colat(i, 0));
        for (int j = 0; j < n_lon; j++) {
            mean_y += pow(y(i, j), 2) * sin(colat(i, j)) * dlon * dlat;
        }
    }
    mean_y /= (4 * M_PI);
    mean_y = sqrt(mean_y);
    ASSERT_NEAR((mean_y - collinear.synthesis_average(functional)) / mean_y, 0,
                1e-5);
}
