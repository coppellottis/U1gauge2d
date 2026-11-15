#ifndef OBSERVABLES_H
    #define OBSERVABLES_H 1

    #include "lattice.h"

    double wilson_loop(Lattice* lat, int R, int T);
    double wilson_loop_mh(Lattice* lat, int R, int T, int N);
    double plaq_angle(Lattice* lat, int s);
    double normalize_angle(double theta);

#endif