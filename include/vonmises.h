#ifndef VONMISES_H
    #define VONMISES_H 1

    #include "pcg_basic.h"
    #include "lattice.h"

    void vonmises_param(Lattice* lat, int s, int mu, double* arr); 
    double vonmises_sample(Lattice* lat, int s, int mu, pcg32_random_t* rng);
#endif