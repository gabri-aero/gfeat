#ifndef _RADIAL_HPP_
#define _RADIAL_HPP_

#include "../Planet.hpp"
#include "BaseObservation.hpp"

class Radial : public BaseObservation {
    double H_(int l, int m, int k) const override final {
        double beta_km = k - m * static_cast<double>(this->Nd) / this->Nr;
        return mu / pow(this->r, 2) * (-(l + 1) * beta_km + 2 * k) /
               (beta_km * pow(this->n, 2) * (1 - pow(beta_km, 2))) *
               pow(ae / this->r, l) * this->flmp.get_Flmk(l, m, k);
    }

public:
    Radial(int l_max, int Nr, int Nd, double I, double we_0 = 0.0,
           double wo_0 = 0.0)
        : BaseObservation(l_max, Nr, Nd, I, we_0, wo_0) {
        this->delta_psi_km_0 = 0;
        // Need to define a constructor to be able to load H function override
        this->initialize_block_observation_system();
    }
};

#endif //_RADIAL_HPP_