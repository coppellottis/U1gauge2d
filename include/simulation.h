#ifndef SIMULATIO_H
    #define SIMULATION_H 1

    #include <stdbool.h>
    #include <stdint.h>
    #include <math.h>
    #include "pcg_basic.h"

    void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng);
    static metropolis_sweep(Lattice* lat, pcg32_random_t* rng);
    static plaq_angle(Lattice* lat, int s);

#endif