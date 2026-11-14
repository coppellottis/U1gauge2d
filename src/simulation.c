#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "lattice.h"
#include "simulation.h"
#include "pcg_basic.h"

static double plaq_angle(Lattice* lat, int s) {
    int sx = lat->up[2*s + 0]; // s + 0
    int sy = lat->up[2*s + 1]; // s + 1

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

            double S_old = -lat->beta*(cos(plaq_angle(lat,p1))+cos(plaq_angle(lat,p2)));

            double eps = 2.5*sqrt(lat->beta);
            double max_delta = (M_PI > eps) ? eps : M_PI;

            double r = (double)pcg32_random_r(rng)/4294967296.0;
            double delta = (2.0*r-1.0)*max_delta;

            double th_old = lat->th[2*s+mu];
            lat->th[2*s+mu] = th_old + delta;

            double S_new = -lat->beta*(cos(plaq_angle(lat,p1))+cos(plaq_angle(lat,p2)));

            double delta_S = S_new - S_old;

            double r_acc = (double)pcg32_random_r(rng) / 4294967296.0;
            if(delta_S>0.0 && r_acc> exp(-delta_S)) {
                lat->th[2*s+mu] = th_old;
            }
        }
    }
    return;
}

void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng) {

    return;
}