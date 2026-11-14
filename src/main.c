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

int main(void) {

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (uint64_t)time(NULL), 42u);  

    Lattice* lat = init_lattice(16, 2, 1, &rng);
    simulation(lat, 1e3, &rng);

    free_lattice(lat);

    return 0;
}