#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "lattice.h"
#include "pcg_basic.h"
#include "simulation.h"
#include "observables.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int main(int argc, char *argv[]) {

    int L = (argc > 1) ? atoi(argv[1]) : 16;
    double beta = (argc > 2) ? atof(argv[2]) : 2;
    int n_sweep = (argc > 3) ? (int)strtod(argv[3], NULL) : (int)1e5;
    const char *filename = (argc > 4) ? argv[4] : "test";
    int R = (argc > 5) ? atoi(argv[5]) : 1;
    int T = (argc > 6) ? atoi(argv[6]) : 1;

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (uint64_t)time(NULL), 42u);  

    Lattice* lat = init_lattice(L, beta, 1, &rng);
    simulation(lat, n_sweep, &rng, filename, R, T);

    free_lattice(lat);

    return 0;
}