#ifndef _MULTI_OBSERVATION_HPP_
#define _MULTI_OBSERVATION_HPP_

#include "AbstractKiteSystem.hpp"
#include "BaseObservation.hpp"

class MultiObservation : public AbstractKiteSystem {
protected:
    std::vector<std::shared_ptr<BaseObservation>> observations;

public:
    MultiObservation(int l_max, int Nr, int Nd,
                     std::vector<std::shared_ptr<BaseObservation>> observations)
        : AbstractKiteSystem(l_max, Nr, Nd), observations(observations) {
        // Compute sampling frequency
        this->df = this->observations.at(0)->get_df();
    }

    auto get_observations() const { return this->observations; }

    void compute_normal_matrix() override final {
        for (auto &observation : this->observations) {
            observation->compute_normal_matrix();
            for (int m0 = 0; m0 < this->n_blocks; m0++) {
                this->N.at(m0) += observation->get_N_blocks().at(m0);
            }
        };
    }
};

#endif // _MULTI_OBSERVATION_HPP_