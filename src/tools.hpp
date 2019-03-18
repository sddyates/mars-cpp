
#pragma once

#include <iostream>
#include <cmath>
#include <hdf5.h>

#include "algorithms.hpp"
#include "grid.hpp"
#include "parameters.hpp"
#include "enumerators.hpp"

struct timeStructure {
    double t;
    double dt;
    double velocity_max;
    double mach_number;
};

inline void prims_to_cons(double *V, double *U, const Grid &g, const Algorithm &a) {

    double v2;

    //std::cout << "V[" << rho << "] = " << V[rho] << std::endl;
    //std::cout << "U[" << rho << "] = " << U[rho] << std::endl;

    U[rho] = V[rho];
    U[mx1] = V[rho]*V[vx1];
    U[mx2] = V[rho]*V[vx2];
    U[mx3] = V[rho]*V[vx3];

    v2 = V[vx1]*V[vx1] + V[vx2]*V[vx2] + V[vx3]*V[vx3];
    U[eng] = 0.5*V[rho]*v2 + V[prs]/a.gamma_1;

    return;
}

inline void cons_to_prims(double *U, double *V, const Grid &g, const Algorithm &a) {

    double m2, kinE;

    m2 = U[mx1]*U[mx1] + U[mx2]*U[mx2] + U[mx3]*U[mx3];
    kinE = 0.5*m2/U[rho];

    if (U[eng] < 0.0) {
        U[eng] = a.small_pressure/a.gamma_1 + kinE;
    }

    V[rho] = U[rho];
    V[vx1] = U[mx1]/U[rho];
    V[vx1] = U[mx2]/U[rho];
    V[vx1] = U[mx3]/U[rho];
    V[prs] = a.gamma_1*(U[eng] - kinE);

    if (V[prs] < 0.0) {
        V[prs] = a.small_pressure;
    }

    return;
}

inline void time_step(double ****V, Grid &g, const Algorithm &a,
    timeStructure time) {

    double cs, cs_max = 0.0;
    double velocity, velocity_max = 0.0;
    int i, j, k;

    for (k=g.kbeg; k<=g.kend; k++) {
        for (j=g.jbeg; j<=g.jend; j++) {
            for (i=g.ibeg; i<=g.iend; i++) {

                cs = std::sqrt(a.gamma*V[prs][k][j][i]/V[rho][k][j][i]);

                if (cs > cs_max) {
                    cs_max = cs;
                }

                velocity = std::max(
                    fabs(V[vx1][k][j][i]), fabs(V[vx2][k][j][i]));
                velocity = std::max(
                    fabs(velocity), fabs(V[vx3][k][j][i]));

                if (velocity > velocity_max) {
                    velocity_max = velocity;
                }

            }
        }
    }

    double speed_max = velocity_max + cs_max;
    double dt = a.cfl*g.dxi_min/speed_max;
    double mach_number = velocity_max/cs_max;

    if (dt < a.small_dt || isnan(dt)) {
        //std::cout << "Error, dt is nan to small, dt = " <<
        //    dt << std::endl;
    }

    dt = fmin(dt, a.max_dt_increase*time.dt);

    if ((time.t + dt) > a.max_time) {
        dt = a.max_time - time.t;
    }

    time.dt = dt;
    time.velocity_max = velocity_max;
    time.mach_number = mach_number;

    return;
}

inline void dump (double ****V, Grid grid, Parameters parameter, int num) {

    hid_t file, dataset, dataspace;
    hsize_t dimsf[4];
    herr_t status;

    file = H5Fcreate("example.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dimsf[0] = grid.nx1_tot;
    dimsf[1] = grid.nx2_tot;
    dimsf[2] = grid.nx3_tot;
    dimsf[3] = grid.nvar;

    //Dataspace creation
    dataspace = H5Screate_simple(2, dimsf, NULL);

    //Dataset creation
    dataset = H5Dcreate(file, "conservative_data", H5T_NATIVE_DOUBLE, dataspace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Actual data IO
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, V);

    //Closing all opened HDF5 objects
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);

    return;
}
