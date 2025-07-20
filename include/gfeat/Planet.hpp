#ifndef _PLANET_HPP_
#define _PLANET_HPP_

struct Planet {
    // Define default planet properties
    double mu = 0.39860044150E+15;
    double theta_dot = 7.2921159e-5;
    double C20 = 1.08263e-3;
    double ae = 0.63781363000E+07;
    double rho_e = 5517;
    double rho_w = 1000;
    double K = 1e-5; // Kaula constant
    // Copy operator
    Planet &operator=(const Planet &other) {
        if (this != &other) {
            this->ae = other.ae;
            // maybe other stuff
        }
        return *this;
    }
};

Planet planet;

// Define alias to planet variables
auto &ae = planet.ae;
auto &mu = planet.mu;
auto &theta_dot = planet.theta_dot;
auto &C20 = planet.C20;
auto &rho_e = planet.rho_e;
auto &rho_w = planet.rho_w;
auto &K = planet.K;

#endif //_PLANET_HPP_
