
#pragma once

#include <algorithm>
#include <math.h>

#include "algorithms.hpp"
#include "enumerators.hpp"
#include "data_column.hpp"

inline void hll(dataColumn &c, const Grid &g, const Algorithm &a, int vxn, int vxt, int vxb)
   /*
    * Estimate the leftmost and rightmost wave signal
    * speeds bounding the Riemann fan based on the
    * input states VL and VR accourding to the Davis
    * Method.
    */
{
    int i, var;
    double csL[c.nx_tot];
    double csR[c.nx_tot];
    double sL_min[c.nx_tot];
    double sL_max[c.nx_tot];
    double sR_min[c.nx_tot];
    double sR_max[c.nx_tot];
    double scrh[c.nx_tot];

    for (i = 0; i < c.nx_tot; i++) {
        csL[i] = sqrt(a.gamma*c.VL[prs][i]/c.VL[rho][i]);
        sL_min[i] = c.VL[vxn][i] - csL[i];
        sL_max[i] = c.VL[vxn][i] + csL[i];

        csR[i] = sqrt(a.gamma*c.VR[prs][i]/c.VR[rho][i]);
        sR_min[i] = c.VR[vxn][i] - csR[i];
        sR_max[i] = c.VR[vxn][i] + csR[i];

        c.SL[i] = std::min(sL_min[i], sR_min[i]);
        c.SR[i] = std::max(sL_max[i], sR_max[i]);

        c.cmax[i] = std::max(fabs(c.SL[i]), fabs(c.SR[i]));
    }

    for (var = 0; var < g.nvar; var++) {
        for (i = 0; i < c.nx_tot; i++) {

            if (c.SL[i] > 0.0) {
                c.flux[var][i] = c.FL[var][i];
                c.pres[i] = c.VL[prs][i];
            }
            else if (c.SR[i] < 0.0) {
                c.flux[var][i] = c.FR[var][i];
                c.pres[i] = c.VR[prs][i];
            }

            else{
                scrh[i] = 1.0/(c.SR[i] - c.SL[i]);
                c.flux[var][i] = c.SL[i]*c.SR[i]*(c.UR[var][i]
                    - c.UL[var][i]) \
                    + c.SR[i]*c.FL[var][i]
                    - c.SL[i]*c.FR[var][i];
                c.flux[var][i] *= scrh[i];
                c.pres[i] = (c.SR[i]*c.VL[prs][i]
                    - c.SL[i]*c.VR[prs][i])*scrh[i];
            }

        }
    }

    return;

}
