#ifndef LATTICE_H
    #define LATTICE_H 1

    #include <stdbool.h>
    #include <stdint.h>
    #include <math.h>
    #include "pcg_basic.h"

    typedef struct{
        int L, V;
        double beta;
        double* th;
        int *up, *dn;
    } Lattice;

    Lattice* init_lattice(int N, double beta, bool hot, pcg32_random_t* rng);
    void free_lattice(Lattice* lat);

#endif