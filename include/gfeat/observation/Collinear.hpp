#ifndef _COLLINEAR_HPP_
#define _COLLINEAR_HPP_

#include "../Planet.hpp"
#include "AlongTrack.hpp"
#include "BaseObservation.hpp"

class Collinear : public BaseObservation {
    double eta;
    double H_(int l, int m, int k) const override final {
        double beta_km = k - m * static_cast<double>(this->Nd) / this->Nr;
        double n = this->n;
        double r = this->r;
        double alpha_lmk =
            ((-(l + 1) * beta_km + 2 * k) * beta_km * cos(beta_km * eta) *
                 sin(eta) +
             (2 * (l + 1) * beta_km - k * (3 + beta_km * beta_km)) *
                 sin(beta_km * eta) * cos(eta)) /
            (n * n * beta_km * beta_km * (1 - beta_km * beta_km));
        return 2 * mu / (r * r) * pow(ae / r, l) * alpha_lmk *
               this->flmp.get_Flmk(l, m, k);
    }

public:
    Collinear(int l_max, double I, double rho_0, int Nr, int Nd,
              double we_0 = 0.0, double wo_0 = 0.0)
        : BaseObservation(l_max, I, Nr, Nd, we_0, wo_0) {
        this->eta = asin(rho_0 / (2 * this->r));
        // Need to define a constructor to be able to load H function override
        this->initialize_block_observation_system();
    }

    using BaseObservation::set_observation_error;

    void
    set_observation_error(std::function<double(double)> range_asd,
                          std::function<double(double)> accelerometer_asd) {
        auto combined_error = [this, range_asd,
                               accelerometer_asd](double f) -> double {
            double w = 2 * M_PI * f;
            double n = this->wo_dot;
            double H_acc_v =
                abs((3 * n * n + w * w) / (w * w * (n * n - w * w)));
            double H_acc_u = abs(2 * n / (w * (n * n - w * w)));
            return range_asd(f) +
                   2 * (H_acc_u + H_acc_v) * accelerometer_asd(f);
        };
        this->set_observation_error(combined_error);
    }
};

#endif //_COLLINEAR_HPP_