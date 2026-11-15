#include <stdlib.h>
#include <math.h>
#include "observables.h"
#include "lattice.h"

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
double wilson_loop_mh(Lattice* lat, int R, int T, int N) {

}