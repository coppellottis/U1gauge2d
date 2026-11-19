#ifndef OBSERVABLES_H
    #define OBSERVABLES_H 1

    #include "lattice.h"
    #include "pcg_basic.h"

    double wilson_loop(Lattice* lat, int R, int T);
    double wilson_loop_mh(Lattice* lat, int R, int T, int N, pcg32_random_t* rng);
    static void multihit_links(Lattice* lat, int mu, int n_hit, pcg32_random_t* rng, double* rho, double* psi);
    double plaq_angle(Lattice* lat, int s);
    double normalize_angle(double theta);
    int topological_charge(Lattice* lat);

#endif