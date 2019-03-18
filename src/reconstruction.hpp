
#pragma once

#include <algorithm>
#include <math.h>

#include "grid.hpp"
#include "data_column.hpp"

inline void minmod(dataColumn &c, double dxi) {
   /*
    * Obtain the left and right hand states from linear extrapolation
    * and apply the minmod limiter.
    */

    int var;
    int i;

    double a;
    double b;
    double dxi_inv = 1.0/dxi;
    double dxi_half = 0.5*dxi;

    double m[c.nx_tot-1+c.gz];

    for (var = 0; var < c.nvar; var++) {
        for (i = c.beg; i < c.end-1; i++) {

            a = (c.V[var][i+1] - c.V[var][i])*dxi_inv;
            b = (c.V[var][i+2] - c.V[var][i+1])*dxi_inv;

            if (a*b < 0.0) {
                m[i] = 0.0;
            }
            else if (fabs(a) < fabs(b)) {
                m[i] = a;
            }
            else{
                m[i] = b;
            }

        }

        for (i = c.beg; i < c.end-1; i++) {
            c.VL[var][i] = c.V[var][i+1] + m[i]*dxi_half;
            c.VR[var][i] = c.V[var][i+2] - m[i+1]*dxi_half;
        }
    }

    return;
}
