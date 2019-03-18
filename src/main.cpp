
/*
 * Main program loop.
 *
 * Simon Daley-Yates @ 01/11/2018
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

#include "parameters.hpp"
#include "grid.hpp"
#include "algorithms.hpp"
#include "tools.hpp"
#include "time_stepping.hpp"
#include "make_arrays.hpp"

int main(int argc, char* argv[]) {
   /*
    * Evolve the (M)HD equations through t = 0 to t = "max time".
    */

    Parameters parameter;
    initilise_parameters(parameter);

    // Set global parameters.
    cout << "" << endl;
    cout << "  -----------------------------------------------" << endl;
    cout << "                                                 " << endl;
    cout << "      //////////      /     |///// ////// ++     " << endl;
    cout << "      ||  ||  ||    ////    ||  // ||            " << endl;
    cout << "      ||  ||  ||   //  //   ||//// //////        " << endl;
    cout << "      ||  ||  ||  ////////  ||  ||     ||        " << endl;
    cout << "      ||  ||  || //      // ||  || ////// 0.2    " << endl;
    cout << "                                                 " << endl;
    cout << "  -----------------------------------------------" << endl;
    cout << "" << endl;
    cout << "  Problem settings:" << endl;
    cout << "      - Name: " << parameter.name << endl;
    cout << "      - Max time: " << parameter.max_time << endl;
    cout << "" << endl;


    // Initialise Algorithms.
    Algorithm algorithm(parameter);
    cout << "  Algorithms initialised..." << endl;

    // Initialise grid.
    Grid grid(parameter);
    cout << "  Grid initialised..." << endl;

    // Array for the primative variables.
    double**** V = MakeArray4D(grid.nvar, grid.nx3_tot,
        grid.nx2_tot, grid.nx1_tot);
    cout << "  State vector created..." << endl;

    // Initialise the state vector accourding to
    // user defined problem.
    // problem.initialise_grid(V, grid);

    int k, j, i;
    for (k = grid.kbeg; k <= grid.kend; k++) {
        for (j = grid.jbeg; j <= grid.jend; j++) {
            for (i = grid.ibeg; i <= grid.iend; i++) {
                if (grid.x1[i] < 0.5) {
                    V[rho][k][j][i] = 1.0;
                    V[prs][k][j][i] = 1.0;
                    V[vx1][k][j][i] = 0.0;
                } else {
                    V[rho][k][j][i] = 0.125;
                    V[prs][k][j][i] = 0.1;
                    V[vx1][k][j][i] = 0.0;
                }
            }
        }
    }
    cout << "    - State vector initialised." << endl;

    // Apply boundary conditions.
    grid.boundary(V, parameter);
    cout << "  Applying boundary conditions..." << endl;

    // First output.
    dump(V, grid, parameter, 0);

    timeStructure time;
    // Perform main integration loop.
    time.t = 0.0;
    time.dt = parameter.initial_dt;
    int n = 0;
    int num = 1;
    double percent = 100.0/parameter.max_time;

    // Integrate in time.
    cout << "  Begining time integration:" << endl;
    cout << "" << endl;
    while (time.t < parameter.max_time) {

        RungaKutta2(V, time.dt, grid, algorithm, parameter);

        time_step(V, grid, algorithm, time);

        cout.precision(2);
        cout << "  n = " << n << "; t = " << time.t << "; dt = "
            << time.dt << "; " << percent*time.t
            << "%; [" << time.velocity_max << ","
            << time.mach_number << "]" << endl;

        if ((time.t + time.dt) > num*parameter.plot_frequency) {
            dump(V, grid, parameter, num);
            num += 1;
        }

        time.t += time.dt;
        n++;

    }

    cout << "  n = " << n << "; t = " << time.t << "; dt = "
        << time.t << "; " << percent*time.t
        << "%; [" << time.velocity_max << ","
        << time.mach_number << "]" << endl;

    RungaKutta2(V, time.dt, grid, algorithm, parameter);
    dump(V, grid, parameter, num);


    cout << "" << endl;
    cout << "  Simulation " << parameter.name
        << " complete..." << endl;
    cout << "" << endl;

}
