#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "lattice.h"
#include "pcg_basic.h"

Lattice* init_lattice(int L, double beta, bool hot, pcg32_random_t* rng) {
    Lattice* lat = malloc(sizeof(Lattice));
    lat->L = L;
    lat->V = L*L;
    lat->beta = beta;

    lat->th = malloc(2*lat->V*sizeof(*lat->th));
    lat->up = malloc(2*lat->V*sizeof(*lat->up));
    lat->dn = malloc(2*lat->V*sizeof(*lat->dn));

    for(int x = 0; x < L; x++) {
        for(int y=0; y<L; y++) {
            int s = x+L*y;

            // Inizializzazione link
            for(int mu = 0; mu < 2; mu++) {
                double rand_01 = (hot)? (double)pcg32_random_r(rng) / 4294967296.0 : 0;
                lat->th[2*s + mu] = (2*rand_01-1)*M_PI;
            }

            // Direzione x
            lat->up[2*s + 0] = (x + 1) % L + L*y;
            lat->dn[2*s + 0] = (x + L - 1) % L + L*y;

            // Direzione y
            lat->up[2*s + 1] = x + L*((y + 1) % L);
            lat->dn[2*s + 1] = x + L*((y + L - 1) % L);
        }
    }
    return lat;
}


void free_lattice(Lattice* lat) {
    free(lat->th);
    free(lat->up);
    free(lat->dn);
    free(lat);
    return;
}
