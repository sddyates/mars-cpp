
#include <iostream>
#include <algorithm>
#include <math.h>

#include "grid.hpp"
#include "make_arrays.hpp"

Grid::Grid (const Parameters &parameter) {

    /*
    if (strcmp(parameter.physics, "hydro") == 0) {
        nvar = 5;
    }
    else {
        nvar = 8;
    }
    */
    nvar = 5;
    gz = 2;

    x1min = parameter.x1_min;
    x1max = parameter.x1_max;
    x2min = parameter.x2_min;
    x2max = parameter.x2_max;
    x3min = parameter.x3_min;
    x3max = parameter.x3_max;

    nx1 = parameter.x1_resolution;
    nx2 = parameter.x2_resolution;
    nx3 = parameter.x3_resolution;

    nx1_tot = nx1 + 2*gz;
    nx2_tot = nx2 + 2*gz;
    nx3_tot = nx3 + 2*gz;

    x1 = MakeArray1D(nx1_tot);
    x2 = MakeArray1D(nx2_tot);
    x3 = MakeArray1D(nx3_tot);

    x1_verts = MakeArray1D(nx1_tot+1);
    x2_verts = MakeArray1D(nx2_tot+1);
    x3_verts = MakeArray1D(nx3_tot+1);

    dx1 = (x1max - x1min)/nx1;
    dx2 = (x2max - x2min)/nx2;
    dx3 = (x3max - x3min)/nx3;
    dxi[0] = dx1; dxi[1] = dx2; dxi[2] = dx3;
    dxi_min = dxi[0];
    for (i = 0; i < 4; i++) {
        if (dxi_min > dxi[i]) {
            dxi_min = dxi[i];
        }
    }

    dv = dx1*dx2*dx3;

    ibeg = gz;
    iend = nx1 + gz;
    jbeg = gz;
    jend = nx2 + gz;
    kbeg = gz;
    kend = nx3 + gz;

    lower_bc_ibeg = 0;
    lower_bc_iend = gz - 1;
    lower_bc_jbeg = 0;
    lower_bc_jend = gz - 1;
    lower_bc_kbeg = 0;
    lower_bc_kend = gz - 1;

    upper_bc_ibeg = nx1 + gz;
    upper_bc_iend = nx1 + 2*gz - 1;
    upper_bc_jbeg = nx2 + gz;
    upper_bc_jend = nx2 + 2*gz - 1;
    upper_bc_kbeg = nx3 + gz;
    upper_bc_kend = nx3 + 2*gz - 1;

    imax = upper_bc_iend;
    jmax = upper_bc_jend;
    kmax = upper_bc_kend;

    cell_centers(x1, x1min, x1max, dx1, nx1_tot);
    cell_centers(x2, x2min, x2max, dx2, nx2_tot);
    cell_centers(x3, x3min, x3max, dx3, nx3_tot);

    cell_edges(x1_verts, x1, dx1, nx1_tot);
    cell_edges(x2_verts, x2, dx2, nx2_tot);
    cell_edges(x3_verts, x3, dx3, nx3_tot);

}

void Grid::boundary(double ****K1, Parameters &p){
    return;
}

void Grid::cell_centers(double* x, double xmin, double xmax,
                        double dx, int nx_tot) {

    int i;
    for (i = 0; i < nx_tot; i++) {
        x[i] = xmin + ((double)i - 1.5)*dx;
    }

    return;
}

void Grid::cell_edges(double* verts, double *x, double dx, int nx_tot) {

    int i;
    for (i = 0; i < nx_tot; i++) {
        verts[i] = x[i] - dx/2.0;
    }
    verts[nx_tot] = x[nx_tot-1] + dx/2.0;

    return;
}
