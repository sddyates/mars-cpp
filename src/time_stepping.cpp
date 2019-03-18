
#include "grid.hpp"
#include "algorithms.hpp"
#include "tools.hpp"
#include "make_arrays.hpp"
#include "right_hand_side.hpp"

void RungaKutta2(double ****V, double dt, Grid &grid,
                    const Algorithm &a, Parameters &p) {

    int k, j, i, var, v2, m2, kinE;

    double**** U = MakeArray4D(grid.nvar, grid.nx3_tot, grid.nx2_tot, grid.nx1_tot);
    double**** K1 = MakeArray4D(grid.nvar, grid.nx3_tot, grid.nx2_tot, grid.nx1_tot);
    double**** K2 = MakeArray4D(grid.nvar, grid.nx3_tot, grid.nx2_tot, grid.nx1_tot);
    double**** U_mid = MakeArray4D(grid.nvar, grid.nx3_tot, grid.nx2_tot, grid.nx1_tot);

    //std::cout << "grid.nvar = " << grid.nvar << std::endl;
    //std::cout << "U = " << U[4][4][11][11] << std::endl;

    //std::cout << "kbeg = " << grid.kbeg << " kend = " << grid.kend << std::endl;
    //std::cout << "jbeg = " << grid.jbeg << " jend = " << grid.jend << std::endl;
    //std::cout << "ibeg = " << grid.ibeg << " iend = " << grid.iend << std::endl;

    for (k = grid.kbeg; k < grid.kend; k++) {
        for (j = grid.jbeg; j < grid.jend; j++) {
            for (i = grid.ibeg; i < grid.iend; i++) {
                //std::cout << "i = " << i << " j = " << j << " k = " << k << std::endl;
                U[rho][k][j][i] = V[rho][k][j][i];
                U[mx1][k][j][i] = V[rho][k][j][i]*V[vx1][k][j][i];
                U[mx2][k][j][i] = V[rho][k][j][i]*V[vx2][k][j][i];
                U[mx3][k][j][i] = V[rho][k][j][i]*V[vx3][k][j][i];

                v2 = V[vx1][k][j][i]*V[vx1][k][j][i]
                   + V[vx2][k][j][i]*V[vx2][k][j][i]
                   + V[vx3][k][j][i]*V[vx3][k][j][i];
                U[eng][k][j][i] = 0.5*V[rho][k][j][i]*v2 + V[prs][k][j][i]/a.gamma_1;
                //prims_to_cons(V, U, grid, a);
            }
        }
    }

    K1 = RHSOperator(U, grid, a);
    grid.boundary(K1, p);

    // May need to recalculate the time step here.

    for (k = grid.kbeg; k < grid.kend; k++) {
        for (j = grid.jbeg; j < grid.jend; j++) {
            for (i = grid.ibeg; i < grid.iend; i++) {
                for (var = 0; var < grid.nvar; var++) {
                    U_mid[var][k][j][i] = U[var][k][j][i]
                                        + dt*K1[var][k][j][i];
                }
            }
        }
    }

    K2 = RHSOperator(U_mid, grid, a);
    grid.boundary(U, p);

    for (k = grid.kbeg; k < grid.kend; k++) {
        for (j = grid.jbeg; j < grid.jend; j++) {
            for (i = grid.ibeg; i < grid.iend; i++) {
                for (var = 0; var < grid.nvar; var++) {
                    U[var][k][j][i] += 0.5*(K1[var][k][j][i]
                                          + dt*K2[var][k][j][i]);
                }
            }
        }
    }

    grid.boundary(U, p);

    for (k = grid.kbeg; k < grid.kend; k++) {
        for (j = grid.jbeg; j < grid.jend; j++) {
            for (i = grid.ibeg; i < grid.iend; i++) {
                //cons_to_prims(U[k][j][i], V[k][j][i], grid, a);

                m2 = U[mx1][k][j][i]*U[mx1][k][j][i]
                   + U[mx2][k][j][i]*U[mx2][k][j][i]
                   + U[mx3][k][j][i]*U[mx3][k][j][i];
                kinE = 0.5*m2/U[rho][k][j][i];

                if (U[eng][k][j][i] < 0.0) {
                    U[eng][k][j][i] = a.small_pressure/a.gamma_1 + kinE;
                }

                V[rho][k][j][i] = U[rho][k][j][i];
                V[vx1][k][j][i] = U[mx1][k][j][i]/U[rho][k][j][i];
                V[vx1][k][j][i] = U[mx2][k][j][i]/U[rho][k][j][i];
                V[vx1][k][j][i] = U[mx3][k][j][i]/U[rho][k][j][i];
                V[prs][k][j][i] = a.gamma_1*(U[eng][k][j][i] - kinE);

                if (V[prs][k][j][i] < 0.0) {
                    V[prs][k][j][i] = a.small_pressure;
                }
            }
        }
    }

    return;

}
