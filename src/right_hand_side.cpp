
#include "algorithms.hpp"
#include "grid.hpp"
#include "tools.hpp"
#include "reconstruction.hpp"
#include "riemann_solvers.hpp"
#include "make_arrays.hpp"

void flux_tensor(dataColumn &c, const Algorithm &a, int vxn, int vxt, int vxb) {
   /*
    * construct the flux tensor from the
    * conservative and primative vaiables.
    */

    int i;

    // could make all the flux columns into vectors so that they know their
    // own lengths.

    for (i = 0; i < c.nx_tot; i++) {
        c.FL[rho][i] = c.UL[rho][i]*c.VL[vxn][i];
        c.FL[vxn][i] = c.UL[rho][i]*c.VL[vxn][i]*c.VL[vxn][i];
        c.FL[vxt][i] = c.UL[rho][i]*c.VL[vxn][i]*c.VL[vxt][i];
        c.FL[vxb][i] = c.UL[rho][i]*c.VL[vxn][i]*c.VL[vxb][i];
        c.FL[eng][i] = c.VL[vxn][i]*(c.UL[eng][i] + c.VL[prs][i]);

        c.FR[rho][i] = c.UR[rho][i]*c.VR[vxn][i];
        c.FR[vxn][i] = c.UR[rho][i]*c.VR[vxn][i]*c.VR[vxn][i];
        c.FR[vxt][i] = c.UR[rho][i]*c.VR[vxn][i]*c.VR[vxt][i];
        c.FR[vxb][i] = c.UR[rho][i]*c.VR[vxn][i]*c.VR[vxb][i];
        c.FR[eng][i] = c.VR[vxn][i]*(c.UR[eng][i] + c.VR[prs][i]);
    }

    return;
}


void face_flux(dataColumn &c, const Grid& g, const Algorithm& a,
                int vxn, int vxt, int vxb)
/*
 * Construct the fluxes through the cell
 * faces normal to the direction of "axis".
 */
{
    int i;
    double v2, m2, kinE;

    for (i = 0; i < c.nx_tot; i++) {
        //cons_to_prims(c.U[i], c.V[i], g, a);

        m2 = c.U[mx1][i]*c.U[mx1][i]
           + c.U[mx2][i]*c.U[mx2][i]
           + c.U[mx3][i]*c.U[mx3][i];
        kinE = 0.5*m2/c.U[rho][i];

        if (c.U[eng][i] < 0.0) {
            c.U[eng][i] = a.small_pressure/a.gamma_1 + kinE;
        }

        c.V[rho][i] = c.U[rho][i];
        c.V[vx1][i] = c.U[mx1][i]/c.U[rho][i];
        c.V[vx1][i] = c.U[mx2][i]/c.U[rho][i];
        c.V[vx1][i] = c.U[mx3][i]/c.U[rho][i];
        c.V[prs][i] = a.gamma_1*(c.U[eng][i] - kinE);

        if (c.V[prs][i] < 0.0) {
            c.V[prs][i] = a.small_pressure;
        }

    }

    minmod(c, g.dxi[vxn-2]);

    for (i = 0; i < c.nx_tot; i++) {
        //prims_to_cons(c.VL[i], c.UL[i], g, a);
        c.UL[rho][i] = c.VL[rho][i];
        c.UL[mx1][i] = c.VL[rho][i]*c.VL[vx1][i];
        c.UL[mx2][i] = c.VL[rho][i]*c.VL[vx2][i];
        c.UL[mx3][i] = c.VL[rho][i]*c.VL[vx3][i];

        v2 = c.VL[vx1][i]*c.VL[vx1][i]
           + c.VL[vx2][i]*c.VL[vx2][i]
           + c.VL[vx3][i]*c.VL[vx3][i];
        c.UL[eng][i] = 0.5*c.VL[rho][i]*v2 + c.VL[prs][i]/a.gamma_1;
        //std::cout << "beg = " << c.beg << " end = " << c.end << " i = " << i << " U = " << c.UR[rho][i] << std::endl;

        //prims_to_cons(c.VR[i], c.UR[i], g, a);
        c.UR[rho][i] = c.VR[rho][i];
        c.UR[mx1][i] = c.VR[rho][i]*c.VR[vx1][i];
        c.UR[mx2][i] = c.VR[rho][i]*c.VR[vx2][i];
        c.UR[mx3][i] = c.VR[rho][i]*c.VR[vx3][i];

        v2 = c.VR[vx1][i]*c.VR[vx1][i]
           + c.VR[vx2][i]*c.VR[vx2][i]
           + c.VR[vx3][i]*c.VR[vx3][i];
        c.UR[eng][i] = 0.5*c.VR[rho][i]*v2 + c.VR[prs][i]/a.gamma_1;

    }

    flux_tensor(c, a, vxn, vxt, vxb);
    flux_tensor(c, a, vxn, vxt, vxb);

    hll(c, g, a, vxn, vxt, vxb);

    return;

}


double**** RHSOperator(double ****U, const Grid& g, const Algorithm& a) {
/*
 * Determine the right hand side of the HD equations.
 */

    int var, k, j, i;
    int vxn, vxt, vxb;

    double**** dflux = MakeArray4D(g.nvar, g.nx3_tot, g.nx2_tot, g.nx1_tot);
    double**** dflux_x1 = MakeArray4D(g.nvar, g.nx3_tot, g.nx2_tot, g.nx1_tot);
    double**** dflux_x2 = MakeArray4D(g.nvar, g.nx3_tot, g.nx2_tot, g.nx1_tot);
    double**** dflux_x3 = MakeArray4D(g.nvar, g.nx3_tot, g.nx2_tot, g.nx1_tot);

    // Loop over column in x1 direction.
    dataColumn c1(g.gz, g.nvar, g.nx1_tot);

    for (k = g.kbeg; k < g.kend; k++) {
        for (j = g.jbeg; j < g.jend; j++) {

            for (var = 0; var < g.nvar; var++) {
                // Need to send the hole column of the normal direction.
                for (i = 0; i < g.nx1_tot; i++) {
                    c1.U[var][i] = U[var][k][j][i];
                }
            }

            face_flux(c1, g, a, vxn=2, vxt=3, vxb=4);

            for (var = 0; var < g.nvar; var++) {
                for (i = g.ibeg; i < g.iend; i++) {
                    //std::cout << "k = " << k << " j = " << j << " var = " << var << " i = " << i << std::endl;
                    dflux_x1[var][k][j][i] = -(c1.flux[var][i+1]
                        - c1.flux[var][i]);
                    dflux_x1[mx1][k][j][i] -= c1.pres[i+1] - c1.pres[i];
                }
            }

        }
    }

    // Loop over column in x2 direction.
    dataColumn c2(g.gz, g.nvar, g.nx2_tot);

    for (k = g.kbeg; k < g.kend; k++) {
        for (i = g.ibeg; i < g.iend; i++) {

            for (var = 0; var < g.nvar; var++) {
                for (j = 0; j < g.nx2_tot; j++) {
                    c2.U[var][j] = U[var][k][j][i];
                }
            }

            face_flux(c2, g, a, vxn=3, vxt=2, vxb=4);

            for (var = 0; var < g.nvar; var++) {
                for (j = g.jbeg; j < g.jend; j++) {
                    dflux_x2[var][k][j][i] = -(c2.flux[var][j+1]
                        - c2.flux[var][j]);
                    dflux_x2[mx2][k][j][i] -= c2.pres[j+1] - c2.pres[j];
                }
            }

        }
    }

    // Loop over column in x3 direction.
    dataColumn c3(g.gz, g.nvar, g.nx3_tot);

    for (j = g.jbeg; j < g.jend; j++) {
        for (i = g.ibeg; i < g.iend; i++) {

            for (var = 0; var < g.nvar; var++) {
                for (k = 0; k < g.nx3_tot; k++) {
                    c3.U[var][k] = U[var][k][j][i];
                }
            }

            face_flux(c3, g, a, vxn=4, vxt=2, vxb=3);

            for (var = 0; var < g.nvar; var++) {
                for (k = g.kbeg; k < g.kend; k++) {
                    dflux_x3[var][k][j][i] = -(c3.flux[var][k-1]
                        - c3.flux[var][k]);
                    dflux_x3[mx3][k][j][i] -= c3.pres[k-1] - c3.pres[k];
                }
            }

        }
    }

    for (j = g.jbeg; j < g.jend; j++) {
        for (i = g.ibeg; i < g.iend; i++) {
            for (k = g.kbeg; k < g.kend; k++) {
                for (var = 0; var < g.nvar; var++) {
                    dflux[var][k][j][i] = dflux_x1[var][k][j][i]/g.dx1
                                        + dflux_x2[var][k][j][i]/g.dx2
                                        + dflux_x3[var][k][j][i]/g.dx3;
                }
            }
        }
    }

    return dflux;
}
