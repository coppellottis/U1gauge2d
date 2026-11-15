#include <stdlib.h>
#include <math.h>
#include "vonmises.h"
#include "lattice.h"
#include "simulation.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static void vonmises_param(Lattice* lat, int s, int mu, double* arr) {
    int p1 = s;
    int p2 = (mu > 0) ? lat->dn[2*s+0] : lat->dn[2*s+1];

    double a1 = plaq_angle(lat, p1);
    double a2 = plaq_angle(lat, p2);

    double sign = (mu == 0) ? -1.0 : 1.0;

    double phi = lat->th[2*s+mu] + sign * 0.5 * (a1 - a2);
    double k   = 2.0 * lat->beta * cos(0.5 * (a1 + a2));

    /* Moriarty richiede k >= 0:  exp(k cos(t-phi)) = exp(|k| cos(t-phi-pi)) se k<0 */
    if (k < 0.0) { k = -k; phi += M_PI; }

    arr[0] = phi;
    arr[1] = k;

    return;
}

/*cfr. http://journals.aps.org/prd/pdf/10.1103/PhysRevD.25.2185*/
double vonmises_sample(Lattice* lat, int s, int mu, pcg32_random_t* rng) {
    double param[2];
    vonmises_param(lat, s, mu, param);

    double phi = param[0];
    double k = param[1];

    double t;

    /**/

    if (k < 1e-12) {
        double r0 = (double)pcg32_random_r(rng) / 4294967296.0;
        t = M_PI * r0;    /* limite k->0: uniforme */
    } else {
        for (;;) {
            double r1 = (double)pcg32_random_r(rng) / 4294967296.0;
            double r2 = (double)pcg32_random_r(rng) / 4294967296.0;

            double x = 1.0 + log1p((1.0-r1)*expm1(-2.0*k))/k;

            if (k*(sin(0.5*M_PI*x)-x- 0.2105137) > log(r2)) {
                t = 0.5 * M_PI * (1.0 - x);
                break;
            }
        }
    }
    double r3 = (double)pcg32_random_r(rng) / 4294967296.0;
    if (r3 < 0.5) t = -t;

    return normalize_angle(phi + t);
}