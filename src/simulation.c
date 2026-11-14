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

static double normalize_angle(double theta) {
    theta = fmod(theta + M_PI, 2*M_PI);
    theta = (theta < 0) ? theta + 2*M_PI : theta;
    return theta - M_PI;
}

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

static void overrelaxation_update(Lattice* lat, int s, int mu) {
    int p1 = s;
    int p2 = (mu > 0) ? lat->dn[2*s+0] : lat->dn[2*s+1];

    double a1 = plaq_angle(lat,p1);
    double a2 = plaq_angle(lat,p2);

    double S_old = cos(a1)+cos(a2);

    double sign = (mu == 0) ? -1.0 : 1.0;
    double theta_bar = sign*(0.5*normalize_angle(a1-a2))+lat->th[2*s+mu];
    double th_old = lat->th[2*s+mu];
    double th_new = normalize_angle(2*theta_bar - th_old);

    lat->th[2*s+mu] = th_new;

    double b1 = plaq_angle(lat, p1);
    double b2 = plaq_angle(lat, p2);

    double S_new = cos(b1)+cos(b2);

    double delta_S = fabs(S_new - S_old);
    if (delta_S > 1e-12) {
        fprintf(stderr, "\n[FATAL ERROR] Overrelaxation: L'azione locale non si è conservata!\n");
        fprintf(stderr, "  -> Sito: %d, Link mu: %d\n", s, mu);
        fprintf(stderr, "  -> Action old:   %.15f\n", S_old);
        fprintf(stderr, "  -> Action new:   %.15f\n", S_new);
        fprintf(stderr, "  -> Differenza:   %e\n", delta_S);
        fprintf(stderr, "Il programma verrà terminato.\n\n");
        
        exit(EXIT_FAILURE);
    }

    return;
} 

static void overrelaxation_sweep(Lattice* lat) {
    for(int s = 0; s<lat->V; s++){
        for(int mu=0; mu<2; mu++) {
            overrelaxation_update(lat, s, mu);
        }
    }
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
            lat->th[2*s+mu] = normalize_angle(th_old + delta);

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


static int hybrid_sweep(Lattice* lat, pcg32_random_t* rng) {
    /* 1 M -> 5 OR */

    int acc = metropolis_sweep(lat, rng);

    for (int k = 0; k < 5; k++) {
        overrelaxation_sweep(lat);
    }

    return acc;
}
void simulation(Lattice* lat, int n_sweep, pcg32_random_t* rng) {

    /* Inizio del timer */
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int acc = 0;
    for(int i = 0; i < n_sweep; i++) {
        acc += hybrid_sweep(lat, rng);
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

    fprintf(output, "L                          = %d\n", lat->L);
    fprintf(output, "V                          = %d\n", lat->V);
    fprintf(output, "N_links                    = %d\n", n_links);
    fprintf(output, "beta                       = %.8f\n", lat->beta);
    fprintf(output, "N_sweeps (1 M + 5 OR)      = %d\n", n_sweep);
    fprintf(output, "N_attempts (M)             = %d\n", n_attempts);
    fprintf(output, "N_accepted (M)             = %d\n", acc);
    fprintf(output, "Acceptance rate (M)        = %.6f\n", acc_rate);
    fprintf(output, "Execution time             = %.6f s\n", elapsed);
    fprintf(output, "========================================\n\n");

    fclose(output);

    return;
}