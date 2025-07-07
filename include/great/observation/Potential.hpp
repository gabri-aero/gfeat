#ifndef _POTENTIAL_HPP_
#define _POTENTIAL_HPP_

#include "../Planet.hpp"
#include "BaseObservation.hpp"

class Potential : public BaseObservation {
    double H_(int l, int m, int k) const override final {
        return mu / this->r * pow(ae / this->r, l) *
               this->flmp.get_Flmk(l, m, k);
    }

public:
    Potential(int l_max, double I, int Nr, int Nd, double we_0 = 0.0,
              double wo_0 = 0.0)
        : BaseObservation(l_max, I, Nr, Nd, we_0, wo_0) {
        this->delta_psi_km_0 = 0;
        // Need to define a constructor to be able to load H function override
        this->initialize_block_observation_system();
    }
};

#endif //_POTENTIAL_HPP_