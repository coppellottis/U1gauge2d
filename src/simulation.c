#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "lattice.h"
#include "simulation.h"
#include "pcg_basic.h"

static double plaq_angle(Lattice* lat, int s) {

    int sx = lat->up[2*s + 0]; // s + x
    int sy = lat->up[2*s + 1]; // s + y

    double theta_p =
          lat->th[2*s  + 0]
        + lat->th[2*sx + 1]
        - lat->th[2*sy + 0]
        - lat->th[2*s  + 1];

    return theta_p;
}

static void metropolis_sweep(Lattice* lat, pcg32_random_t* rng) {
    int V = (lat->L)*(lat->L);
    for(int s=0; s<V; s++){
        for(int mu=0; mu<2; mu++) {
            int p1 = s;
            int p2 = (mu>0) ? lat->dn[2*s+0] : lat->dn[2*s+1];

            S_old =

        }
    }
    return;
}

void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng) {

    return;
}