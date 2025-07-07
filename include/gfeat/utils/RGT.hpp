#ifndef _UTILS_RGT_HPP_
#define _UTILS_RGT_HPP_

double raan_dot(double sma, double ecc, double inc) {
    // Convert to canonical Delaunay elements
    double L = sqrt(mu * sma);
    double G = L / sqrt(1 - ecc * ecc);
    double H = G * cos(inc);
    // Apply second-order Brouwer theory
    double k2 = -ae * ae * C20 / 2;
    double n0 = pow(mu, 2) / pow(L, 3);
    double gamma2 = pow(mu, 2) * k2 / pow(L, 4);
    double dh_dt = n0 * (-3.0 * gamma2 * pow(L, 4) / pow(G, 4) * H / G +
                         pow(gamma2, 2) * ((27.0 / 8 * pow(L, 6) / pow(G, 6) +
                                            9.0 / 2 * pow(L, 7) / pow(G, 7) -
                                            15.0 / 8 * pow(L, 8) / pow(G, 8)) *
                                               H / G +
                                           (-15.0 / 8 * pow(L, 6) / pow(G, 6) -
                                            27.0 / 2 * pow(L, 7) / pow(G, 7) -
                                            105.0 / 8 * pow(L, 8) / pow(G, 8)) *
                                               pow(H, 3) / pow(G, 3)));
    return dh_dt;
}

double ma_dot(double sma, double ecc, double inc) {
    // Convert to canonical Delaunay elements
    double L = sqrt(mu * sma);
    double G = L / sqrt(1 - ecc * ecc);
    double H = G * cos(inc);
    // Apply second-order Brouwer theory
    double k2 = -ae * ae * C20 / 2;
    double n0 = pow(mu, 2) / pow(L, 3);
    double gamma2 = pow(mu, 2) * k2 / pow(L, 4);
    double dl_dt =
        n0 * (1 +
              3 * gamma2 * pow(L, 3) / pow(G, 3) *
                  (-1.0 / 2 + 3.0 / 2 * pow(H, 2) / pow(G, 2)) +
              pow(gamma2, 2) * (75.0 / 32 * pow(L, 5) / pow(G, 5) +
                                3.0 / 2 * pow(L, 6) / pow(G, 6) -
                                -45.0 / 32 * pow(L, 7) / pow(G, 7) +
                                (-135.0 / 16 * pow(L, 5) / pow(G, 5) -
                                 9 * pow(L, 6) / pow(G, 6) +
                                 45.0 / 16 * pow(L, 7) / pow(G, 7)) *
                                    pow(H, 2) / pow(G, 2) +
                                (75.0 / 32 * pow(L, 5) / pow(G, 5) +
                                 27.0 / 2 * pow(L, 6) / pow(G, 6) +
                                 315.0 / 32 * pow(L, 7) / pow(G, 7)) *
                                    pow(H, 4) / pow(G, 4)));
    return dl_dt;
}

double aop_dot(double sma, double ecc, double inc) {
    // Convert to canonical Delaunay elements
    double L = sqrt(mu * sma);
    double G = L / sqrt(1 - ecc * ecc);
    double H = G * cos(inc);
    // Apply second-order Brouwer theory
    double k2 = -ae * ae * C20 / 2;
    double n0 = pow(mu, 2) / pow(L, 3);
    double gamma2 = pow(mu, 2) * k2 / pow(L, 4);
    double dg_dt =
        n0 * (3 * gamma2 * pow(L, 4) / pow(G, 4) *
                  (-1.0 / 2 + 5.0 / 2 * pow(H, 2) / pow(G, 2)) +
              pow(gamma2, 2) * (75.0 / 32 * pow(L, 6) / pow(G, 6) +
                                9.0 / 4 * pow(L, 7) / pow(G, 7) -
                                105.0 / 32 * pow(L, 8) / pow(G, 8) +
                                (-189.0 / 16 * pow(L, 6) / pow(G, 6) -
                                 18 * pow(L, 7) / pow(G, 7) +
                                 135.0 / 16 * pow(L, 8) / pow(G, 8)) *
                                    pow(H, 2) / pow(G, 2) +
                                (135.0 / 32 * pow(L, 6) / pow(G, 6) +
                                 135.0 / 4 * pow(L, 7) / pow(G, 7) +
                                 1155.0 / 32 * pow(L, 8) / pow(G, 8)) *
                                    pow(H, 4) / pow(G, 4)));
    return dg_dt;
}

double sma_from_rgt(int Nr, int Nd, double sma_0, double inc, double ecc = 0,
                    int max_iter = 100, double h = 1.0, double tol = 1e-3) {
    auto f = [Nr, Nd, ecc, inc](double sma) -> double {
        double we_dot = raan_dot(sma, ecc, inc) - theta_dot;
        double wo_dot = aop_dot(sma, ecc, inc) + ma_dot(sma, ecc, inc);
        return static_cast<double>(Nr) / Nd + wo_dot / we_dot;
    };

    for (int iter = 0; iter < max_iter; iter++) {
        double f_next = f(sma_0);
        double df_next = (f(sma_0 + h) - f(sma_0 - h)) / (2 * h);
        sma_0 -= f_next / df_next;
        if (std::abs(f_next / df_next) < tol)
            break;
    }

    return sma_0;
}

#endif // _UTILS_RGT_HPP_