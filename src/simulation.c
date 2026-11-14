#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include "lattice.h"
#include "simulation.h"
#include "pcg_basic.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


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

static int metropolis_sweep(Lattice* lat, pcg32_random_t* rng) {
    int acc = 0;
    int V = lat->V;
    for(int s=0; s<V; s++){
        for(int mu=0; mu<2; mu++) {
            int p1 = s;
            int p2 = (mu>0) ? lat->dn[2*s+0] : lat->dn[2*s+1];

            double S_old = -lat->beta*(cos(plaq_angle(lat,p1))+cos(plaq_angle(lat,p2)));

            double eps = 2.5/sqrt(lat->beta);
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
            } else acc++;
        }
    }
    return acc;
}

void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng) {

    /* Inizio del timer */
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int acc = 0;
    for(int i = 0; i < n_sweep; i++) {
        acc += metropolis_sweep(lat, rng);
    }

    /* Fine del timer */
    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed =
        (double)(end.tv_sec - start.tv_sec)
        + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;

    /* Statistiche */
    int n_links = 2 * lat->V;
    int n_attempts = n_sweep * n_links;

    double acc_rate = (double)acc / (double)n_attempts;

    /* Apertura del file in append */
    FILE* output = fopen("simulation.log", "a");

    if (output == NULL) {
        perror("Errore nell'apertura del file di log");
        return;
    }

    /* Data e ora */
    time_t now = time(NULL);
    struct tm* t = localtime(&now);

    fprintf(output, "========================================\n");
    fprintf(output, "Date            = %04d-%02d-%02d %02d:%02d:%02d\n",
            t->tm_year + 1900,
            t->tm_mon + 1,
            t->tm_mday,
            t->tm_hour,
            t->tm_min,
            t->tm_sec);

    fprintf(output, "L               = %d\n", lat->L);
    fprintf(output, "V               = %d\n", lat->V);
    fprintf(output, "N_links         = %d\n", n_links);
    fprintf(output, "beta            = %.8f\n", lat->beta);
    fprintf(output, "N_sweeps        = %d\n", n_sweep);
    fprintf(output, "N_attempts      = %d\n", n_attempts);
    fprintf(output, "N_accepted      = %d\n", acc);
    fprintf(output, "Acceptance rate = %.6f\n", acc_rate);
    fprintf(output, "Execution time  = %.6f s\n", elapsed);
    fprintf(output, "========================================\n\n");

    fclose(output);

    return;
}