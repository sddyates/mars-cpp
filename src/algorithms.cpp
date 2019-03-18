
#include <iostream>

#include "algorithms.hpp"
#include "parameters.hpp"
#include "reconstruction.hpp"
#include "riemann_solvers.hpp"

Algorithm::Algorithm (const Parameters &parameter) {

    /*
    this->assign_riemann_solver(parameter);
    this->assign_reconstruction(parameter);
    this->assign_time_stepping(parameter);
    */

    if (strcmp(parameter.dimensions, "1D") == 0 &&
        strcmp(parameter.dimensions, "2D") != 0 &&
        strcmp(parameter.dimensions, "3D") != 0) {
        is_1D = true;
        is_2D = false;
        is_3D = false;
    }
    else if (strcmp(parameter.dimensions, "1D") != 0 &&
        strcmp(parameter.dimensions, "2D") == 0 &&
        strcmp(parameter.dimensions, "3D") == 0) {
        is_1D = false;
        is_2D = true;
        is_3D = false;
    }
    else {
        is_1D = false;
        is_2D = false;
        is_3D = true;
    }

    gamma = parameter.gamma;
    gamma_1 = gamma - 1.0;
    cfl = parameter.cfl;
    small_pressure = 1.0e-12;
    max_dt_increase = parameter.max_dt_increase;
    small_dt = 1.0e-15;

    return;
}

/*
void Algorithm::assign_riemann_solver (const Parameters parameter) {

    // This method assigns the function call for
    // the Riemann solver.

    if (parameter.riemann == "tvdlf") {
        riemann_solver = tvdlf;
    }
    else if (parameter.riemann == "hll") {
        riemann_solver = hll;
    }
    else if (parameter.riemann == "hllc") {
        riemann_solver = hllc;
    }
    else {
        std::cout << "Error: invalid riennman solver." << std::endl;
    }

    return;
}

void Algorithm::assign_reconstruction (const Parameters parameter) {

    // This method assigns the function call for
    // the reconstruction stage.

    if (parameter.reconstruction == "flat") {
        reconstruction = flat;
    }
    else if (parameter.reconstruction == "linear") {
        reconstruction = minmod;
    }
    else {
        std::cout << "Error: Invalid reconstructor." << std::endl;
    }

    return;
}

void Algorithm::assign_time_stepping (const Parameters parameter) {

    // This method assigns the function call for
    // the method used for time stepping.

    if (parameter.time stepping == "Euler") {
        time_incriment = Euler;
    }
    else if (parameter.time_stepping == "RK2") {
        time_incriment = RungaKutta2;
    }
    else {
        std::cout << "Error: Invalid integrator." << std::endl;
    }

    return;
}
*/
