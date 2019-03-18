
#include "globals.cpp"
#include "grid.hpp"
#include "problem.hpp"

void Problem::Problem(){
    this->initilise_parameters();
}

void Problem::initilise_parameters() {

    parameter->name = "shock tube";

    parameter->dimensions = "1D";
    parameter->x1_min = 0.0;
    parameter->x2_min = 1.0;
    parameter->x1_resolution = 128;
    parameter->x2_resolution = 128;
    parameter->x3_resolution = 1;

    parameter->cfl = 0.3;
    parameter->initial_dt = 1.0e-5;
    parameter->max_dt_increase = 1.5;
    parameter->max_time = 1.0;

    parameter->plot_frequency = 0.1;
    parameter->print_to_file = false;

    parameter->gamma = 1.666666;
    parameter->unit_density = 1.0;
    parameter->unit_length = 1.0;
    parameter->unit_velocity = 1.0;

    parameter->riemann = "tvdlf";
    parameter->reconstruction = "linear";
    parameter->limiter = "minmod";
    parameter->time_stepping = "RK2";
    parameter->physics = "hd";

    parameter->lower_x1_boundary = "outflow";
    parameter->lower_x2_boundary = "outflow";
    parameter->upper_x1_boundary = "outflow";
    parameter->upper_x2_boundary = "outflow";
    parameter->internal_boundary = false;

    return;
}

void Problem::initialise_grid(double *V, const Grid& g) {
/*
 * Populate state vector with fluid values.
 */

    for (k = g.kbeg; k = g.kend; k++) {
        for (j = g.jbeg; j = g.jend; j++) {
            for (i = g.ibeg; i = g.iend; i++) {
                if (g.x1[i] < 0.5) {
                    V[rho][k][j][i] = 1.0
                    V[prs][k][j][i] = 1.0
                    V[vx1][k][j][i] = 0.0
                } else {
                    V[rho][k][j][i] = 0.125
                    V[prs][k][j][i] = 0.1
                    V[vx1][k][j][i] = 0.0
                }
            }
        }
    }

}
