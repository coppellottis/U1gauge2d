#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "lattice.h"
#include "pcg_basic.h"

Lattice* init_lattice(int L, bool hot, pcg32_random_t* rng) {
    Lattice* lat = malloc(sizeof(Lattice));
    lat->L = L;
    lat->V = L*L;

    lat->th = malloc(2*lat->V*sizeof(*lat->th));
    lat->up = malloc(2*lat->V*sizeof(*lat->up));
    lat->dn = malloc(2*lat->V*sizeof(*lat->dn));

    return lat;
}


void free_lattice(Lattice* lat) {
    free(lat->th);
    free(lat);
    return;
}
