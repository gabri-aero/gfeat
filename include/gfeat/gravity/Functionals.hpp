#ifndef _FUNCTIONALS_HPP_
#define _FUNCTIONALS_HPP_

#include "../Planet.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <map>

class BaseFunctional {
  public:
    double common_factor;
    BaseFunctional(double common_factor) : common_factor(common_factor) {};
    virtual double degree_dependent_factor(int l) const = 0;
};

class GeoidHeight : public BaseFunctional {
  public:
    GeoidHeight() : BaseFunctional(ae) {};
    double degree_dependent_factor(int l) const override final { return 1; }
};

class GravityAnomaly : public BaseFunctional {
  public:
    GravityAnomaly() : BaseFunctional(1e5 * mu / pow(ae, 2)) {};
    double degree_dependent_factor(int l) const override final { return l - 1; }
};

class EquivalentWaterHeight : public BaseFunctional {
    // Love numbers data from Wahr et al. (1998)
    std::map<int, double> k_l = {
        {0, +0.000},  {1, +0.027},  {2, -0.303},  {3, -0.194},   {4, -0.132},
        {5, -0.104},  {6, -0.089},  {7, -0.081},  {8, -0.076},   {9, -0.072},
        {10, -0.069}, {12, -0.064}, {15, -0.058}, {20, -0.051},  {30, -0.040},
        {40, -0.033}, {50, -0.027}, {70, -0.020}, {100, -0.014}, {150, -0.010},
        {200, -0.007}};

    double smoothing_radius;
    Eigen::VectorXd Wl;

    double get_k_l(int l) const {
        // Find equal or greater than l data in the map;
        auto it = k_l.lower_bound(l);
        // Handle exactly equal case
        if (it->first == l)
            return it->second;
        // Handle case that requires interpolation
        auto upper = it;
        auto lower = std::prev(it);
        // Get values
        int l1 = lower->first;
        int l2 = upper->first;
        double k1 = lower->second;
        double k2 = upper->second;
        // Return interpolation result
        return k1 + (k2 - k1) * static_cast<double>(l - l1) / (l2 - l1);
    }

  public:
    EquivalentWaterHeight(double smoothing_radius)
        : BaseFunctional(ae * rho_e / (3 * rho_w)),
          smoothing_radius(smoothing_radius) {
        int L = 200;
        Wl.resize(L + 1);
        if (smoothing_radius == 0) {
            Wl.setConstant(1);
        } else {
            // Compute b term
            double b = log(2) / (1 - cos(smoothing_radius / ae));
            // Continued fraction algorithm (Gautschi, 1981 and Piretzidis,
            // 2019)
            Eigen::VectorXd rho(L + 1);
            Eigen::VectorXd Wl_hat(L + 1);
            rho(L) = 0;
            // Backward step
            for (int l = L; l >= 1; l--) {
                rho(l - 1) = -1 / (-1 / b - 2 * l / b - rho(l));
            }
            // Forward step
            Wl_hat(0) = 1;
            for (int n = 1; n <= L; n++) {
                Wl_hat(n) = rho(n - 1) * Wl_hat(n - 1);
            }
            // No normalization
            Wl = Wl_hat;
        }
    };
    double degree_dependent_factor(int l) const override final {
        return (2 * l + 1) * Wl(l) / (1 + get_k_l(l));
    }
};

#endif // FUNCTIONALS