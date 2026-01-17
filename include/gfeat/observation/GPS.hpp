#ifndef _GPS_HPP_
#define _GPS_HPP_

#include "AlongTrack.hpp"
#include "CrossTrack.hpp"
#include "MultiObservation.hpp"
#include "Radial.hpp"

class GPS : public MultiObservation {
private:
    static auto get_observations(int l_max, double I, int Nr, int Nd,
                                 double we_0, double wo_0) {
        std::vector<std::shared_ptr<BaseObservation>> gps_observations;
        gps_observations.push_back(
            std::make_shared<Radial>(l_max, Nr, Nd, I, we_0, wo_0));
        gps_observations.push_back(
            std::make_shared<AlongTrack>(l_max, Nr, Nd, I, we_0, wo_0));
        gps_observations.push_back(
            std::make_shared<CrossTrack>(l_max, Nr, Nd, I, we_0, wo_0));
        return gps_observations;
    }

public:
    GPS(int l_max, int Nr, int Nd, double I, double we_0 = 0.0,
        double wo_0 = 0.0)
        : MultiObservation(l_max, Nr, Nd,
                           get_observations(l_max, I, Nr, Nd, we_0, wo_0)) {}

    void set_observation_error(double sigma_u, double sigma_v, double sigma_w,
                               std::function<double(double)> ddu_asd,
                               std::function<double(double)> ddv_asd,
                               std::function<double(double)> ddw_asd) {
        // Radial
        auto radial = this->observations.at(0);
        radial->set_observation_error(
            [radial, ddu_asd, ddv_asd, sigma_u](double f) -> double {
                double w = 2 * M_PI * f;
                double n = radial->get_wo_dot();
                double H_ddu_u = abs(w / (w * (n * n - w * w)));
                double H_ddv_u = abs(2 * n / (w * (n * n - w * w)));
                return sqrt(pow(sigma_u, 2) + pow(H_ddu_u * ddu_asd(f), 2) +
                            pow(H_ddv_u * ddv_asd(f), 2));
            });
        // Along-track
        auto along = this->observations.at(1);
        along->set_observation_error(
            [along, ddu_asd, ddv_asd, sigma_v](double f) -> double {
                double w = 2 * M_PI * f;
                double n = along->get_wo_dot();
                double H_ddu_v = abs(2 * n * w / (w * w * (n * n - w * w)));
                double H_ddv_v =
                    abs((w * w + 3 * n * n) / (w * w * (n * n - w * w)));
                return sqrt(pow(sigma_v, 2) + pow(H_ddu_v * ddu_asd(f), 2) +
                            pow(H_ddv_v * ddv_asd(f), 2));
            });
        // Cross-track
        auto cross = this->observations.at(2);
        cross->set_observation_error(
            [cross, ddw_asd, sigma_w](double f) -> double {
                double w = 2 * M_PI * f;
                double n = cross->get_wo_dot();
                double H_ddw_w = abs(1 / (n * n - w * w));
                return sqrt(pow(sigma_w, 2) + pow(H_ddw_w * ddw_asd(f), 2));
            });
    }
    void set_observation_error(std::function<double(double)> asd) {
        for (auto &observation : this->observations) {
            observation->set_observation_error(asd);
        }
    }
};

#endif // _GPS_HPP_