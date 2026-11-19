#include <stdlib.h>
#include <math.h>
#include "observables.h"
#include "lattice.h"
#include "vonmises.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double normalize_angle(double theta) {
    theta = fmod(theta + M_PI, 2*M_PI);
    theta = (theta < 0) ? theta + 2*M_PI : theta;
    return theta - M_PI;
}

double plaq_angle(Lattice* lat, int s) {
    int sx = lat->up[2*s + 0]; // s + 0
    int sy = lat->up[2*s + 1]; // s + 1

    double theta_p =
          lat->th[2*s  + 0]
        + lat->th[2*sx + 1]
        - lat->th[2*sy + 0]
        - lat->th[2*s  + 1];

    return theta_p;
}

int topological_charge(Lattice* lat) {
    double sum_phi = 0.0;

    for (int s = 0; s < lat->V; s++) {
        /* plaq_angle restituisce l'angolo NON ridotto: va riportato in [-pi,pi) */
        sum_phi += normalize_angle(plaq_angle(lat, s));
    }

    double Q = sum_phi / (2.0 * M_PI);
    int Q_int = (int)lround(Q);

    #ifdef DEBUG
        if (fabs(Q - Q_int) > 1e-8) {
            fprintf(stderr, "[FATAL ERROR] Q non intera: Q = %.12f\n", Q);
            exit(EXIT_FAILURE);
        }
    #endif

    return Q_int;
}

double wilson_loop(Lattice* lat, int R, int T) {
    double sum_cos_th = 0.0;
    int V = lat->V;

    for(int s = 0; s<V; s++) {
        double sum_th = 0.0;
        int curr = s;

        for(int r = 0; r<R; r++) {
            sum_th += lat->th[2*curr+0];
            curr = lat->up[2*curr+0];
        }

        for(int t=0; t<T; t++) {
            sum_th += lat->th[2*curr+1];
            curr = lat->up[2*curr+1];
        }

        for(int r=0; r<R; r++) {
            curr = lat->dn[2*curr+0];
            sum_th -= lat->th[2*curr];
        }

        for(int t=0; t<T; t++) {
            curr = lat->dn[2*curr+1];
            sum_th -= lat->th[2*curr+1];
        }

        sum_cos_th += cos(sum_th);
    }

    return sum_cos_th/V;
}

/* Multi-hit*/

static void multihit_links(Lattice* lat, int mu, int n_hit, pcg32_random_t* rng, double* rho, double* psi) {
    double param[2];

    for(int s = 0; s<lat->V; s++) {
        // vonmises_param(lat, s, mu, param); 

        // double phi = param[0];
        // double k = param[1];

        double X = 0.0, Y = 0.0;
        for(int i = 0; i<n_hit; i++) {
            double th = vonmises_sample(lat, s, mu, rng);

            X += cos(th);
            Y += sin(th);
        }
        
        rho[s] = hypot(X, Y) / n_hit; 
        psi[s] = atan2(Y, X);
    }

    return;
}

double wilson_loop_mh(Lattice* lat, int R, int T, int n_hit, pcg32_random_t* rng) {
    double sum_cos_th = 0.0;
    int V = lat->V;

    double* rho = malloc(V*sizeof(double));
    double* psi = malloc(V*sizeof(double));

    multihit_links(lat, 1, n_hit, rng, rho, psi);

    for (int s = 0; s < V; s++) {
        double sum_th = 0.0;
        double prod   = 1.0;                    
        int curr = s;

        for (int r = 0; r < R; r++) {           /* lato spaziale: invariato */
            sum_th += lat->th[2*curr+0];
            curr = lat->up[2*curr+0];
        }
        for (int t = 0; t < T; t++) {           /* lato temporale in avanti */
            prod   *= rho[curr];
            sum_th += psi[curr];
            curr = lat->up[2*curr+1];
        }
        for (int r = 0; r < R; r++) {           /* lato spaziale: invariato */
            curr = lat->dn[2*curr+0];
            sum_th -= lat->th[2*curr+0];
        }
        for (int t = 0; t < T; t++) {           /* lato temporale all'indietro */
            curr = lat->dn[2*curr+1];
            prod   *= rho[curr];
            sum_th -= psi[curr];
        }

        sum_cos_th += prod * cos(sum_th);
    }

    free(rho);
    free(psi);

    return sum_cos_th / V;
}