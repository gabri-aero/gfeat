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
};

#endif // _GPS_HPP_