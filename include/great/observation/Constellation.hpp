#ifndef _CONSTELLATION_HPP_
#define _CONSTELLATION_HPP_

#include "BaseMultiObservation.hpp"
#include "Collinear.hpp"

enum class LongitudePolicy { INTERLEAVING = 0, OVERLAPPING = 1 };

class Constellation : public MultiObservation {
private:
    static auto get_observations(int l_max, Eigen::VectorXd inclinations,
                                 double rho, int Nr, int Nd,
                                 LongitudePolicy policy) {
        std::vector<std::shared_ptr<BaseObservation>> collinear_observations;
        int Np = inclinations.size();
        int i = (Nr * Nd) % 2;
        int j = policy == LongitudePolicy::INTERLEAVING ? 0 : Np - 1;
        for (int k = 0; k < Np; k++) {
            double dlam_0 = M_PI * (i + 1) * (1 + j) * k / (Nr * Np);
            collinear_observations.push_back(std::make_shared<Collinear>(
                l_max, inclinations[k], rho, Nr, Nd, dlam_0, 0));
        }
        return collinear_observations;
    }

public:
    Constellation(int l_max, Eigen::VectorXd inclinations, double rho, int Nr,
                  int Nd,
                  LongitudePolicy policy = LongitudePolicy::INTERLEAVING)
        : MultiObservation(
              l_max, Nr, Nd,
              get_observations(l_max, inclinations, rho, Nr, Nd, policy)) {}

    void
    set_observation_error(std::function<double(double)> range_asd,
                          std::function<double(double)> accelerometer_asd) {
        for (auto &observation : this->observations) {
            double n = observation->get_wo_dot();
            auto combined_error = [n, range_asd,
                                   accelerometer_asd](double f) -> double {
                double w = 2 * M_PI * f;
                double H_acc_v =
                    abs((3 * n * n + w * w) / (w * w * (n * n - w * w)));
                double H_acc_u = abs(2 * n / (w * (n * n - w * w)));
                return range_asd(f) +
                       2 * (H_acc_u + H_acc_v) * accelerometer_asd(f);
            };
            observation->set_observation_error(combined_error);
        }
        this->update();
    }
};

class Bender : public Constellation {
public:
    Bender(int l_max, Eigen::VectorXd inclinations, double eta, int Nr, int Nd,
           LongitudePolicy policy = LongitudePolicy::INTERLEAVING)
        : Constellation(l_max, inclinations, eta, Nr, Nd, policy) {}
};

#endif // _CONSTELLATION_HPP_