#ifndef _ALONG_TRACK_HPP_
#define _ALONG_TRACK_HPP_

#include "../Planet.hpp"
#include "BaseObservation.hpp"

class AlongTrack : public BaseObservation {
    double H_(int l, int m, int k) const override final {
        double beta_km = k - m * static_cast<double>(this->Nd) / this->Nr;
        return mu / pow(this->r, 2) *
               (-(l + 1) * 2 * beta_km + (3 + pow(beta_km, 2)) * k) /
               (pow(beta_km * this->n, 2) * (1 - pow(beta_km, 2))) *
               pow(ae / this->r, l) * this->flmp.get_Flmk(l, m, k);
    }

public:
    AlongTrack(int l_max, int Nr, int Nd, double I, double we_0 = 0.0,
               double wo_0 = 0.0)
        : BaseObservation(l_max, Nr, Nd, I, we_0, wo_0) {
        this->delta_psi_km_0 = M_PI / 2; // (Akm, Bkm) -> (Slm, -Clm)
        // Need to define a constructor to be able to load H function override
        this->initialize_block_observation_system();
    }
};

#endif //_ALONG_TRACK_HPP_