#ifndef SIMULATION_H
    #define SIMULATION_H 1

    #include <stdbool.h>
    #include <stdint.h>
    #include <math.h>
    #include "pcg_basic.h"

    void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng, const char* filename, int R, int T);
    static int metropolis_sweep(Lattice* lat, pcg32_random_t* rng);
    static int hybrid_sweep(Lattice* lat, pcg32_random_t* rng);
    static void overrelaxation_sweep(Lattice* lat);
    static void overrelaxation_update(Lattice* lat, int s, int mu);
    static void heat_bath_sweep(Lattice* lat, pcg32_random_t* rng);
#endif