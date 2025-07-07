#ifndef _BASE_MULTI_OBSERVATION_HPP_
#define _BASE_MULTI_OBSERVATION_HPP_

#include "AbstractKiteSystem.hpp"
#include "BaseObservation.hpp"

class MultiObservation : public AbstractKiteSystem {
protected:
    std::vector<std::shared_ptr<BaseObservation>> observations;

    Eigen::MatrixXd get_H() override final {
        return observations.at(0)->get_H();
    };
    Eigen::MatrixXd get_Pyy() override final {
        return observations.at(0)->get_Pyy();
    };

public:
    MultiObservation(int l_max, int Nr, int Nd,
                     std::vector<std::shared_ptr<BaseObservation>> observations)
        : AbstractKiteSystem(l_max, Nr, Nd), observations(observations) {
        // Compute sampling frequency
        this->df = this->observations.at(0)->get_df();
    }

    void update() {
        // First define observation sizes
        for (int m0 = 0; m0 < this->n_blocks; m0++) {
            // Compute number of observations per block
            this->ny[m0] = 0; // Initialize counter
            for (auto observation : this->observations) {
                this->ny[m0] += observation->get_ny()[m0];
            }
            // Pre-allocate matrices accordingly
            int nx = this->nx[m0];
            this->H[m0].resize(this->ny[m0], nx);
            this->Pyy[m0].resize(this->ny[m0]);
            // Loop over observations to assign matrices
            int i = 0; // Offset for block
            for (auto observation : this->observations) {
                this->H[m0].block(i, 0, observation->get_ny()[m0], nx) =
                    observation->get_H_blocks()[m0];
                this->Pyy[m0].diagonal().segment(i, observation->get_ny()[m0]) =
                    observation->get_Pyy_blocks()[m0].diagonal();
                i += observation->get_ny()[m0];
            }
        }
    }

    auto get_observations() const { return this->observations; }
};

#endif // _BASE_MULTI_OBSERVATION_HPP_