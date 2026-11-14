#ifndef SIMULATIO_H
    #define SIMULATION_H 1

    #include <stdbool.h>
    #include <stdint.h>
    #include <math.h>
    #include "pcg_basic.h"

    void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng);
    static int metropolis_sweep(Lattice* lat, pcg32_random_t* rng);
    static double plaq_angle(Lattice* lat, int s);
    static double normalize_angle(double theta);
    static int hybrid_sweep(Lattice* lat, pcg32_random_t* rng);
    static void overrelaxation_sweep(Lattice* lat);
    static void overrelaxation_update(Lattice* lat, int s, int mu);
#endif