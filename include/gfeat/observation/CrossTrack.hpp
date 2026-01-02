#ifndef _CROSS_TRACK_HPP_
#define _CROSS_TRACK_HPP_

#include "../Planet.hpp"
#include "BaseObservation.hpp"

class CrossTrack : public BaseObservation {
    double H_(int l, int m, int k) const override final {
        double beta_km = k - m * static_cast<double>(this->Nd) / this->Nr;
        return mu / pow(this->r, 2) * 1 /
               (pow(this->n, 2) * (1 - pow(beta_km, 2))) *
               pow(ae / this->r, l) * this->flmp.get_Flmk_star(l, m, k);
    }

public:
    // Need to re-define starting index
    int l_min(int m, int k) const override final {
        return std::min(BaseObservation::l_min(m, k - 1),
                        BaseObservation::l_min(m, k + 1));
    }

    void compute_flmp() override final {
        this->flmp = Flmp(this->l_max, this->I, true);
    }

    CrossTrack(int l_max, int Nr, int Nd, double I, double we_0 = 0.0,
               double wo_0 = 0.0)
        : BaseObservation(l_max, Nr, Nd, I, we_0, wo_0) {
        this->delta_psi_km_0 = M_PI / 2; // (Akm, Bkm) -> (Slm, -Clm)
        // Need to define a constructor to be able to load H function override
        this->initialize_block_observation_system();
    }
};

#endif //_CROSS_TRACK_HPP_