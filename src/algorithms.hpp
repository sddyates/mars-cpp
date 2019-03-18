
#pragma once

#include "parameters.hpp"

class Algorithm {
public:

    Algorithm(const Parameters &parameter);

    bool is_1D;
    bool is_2D;
    bool is_3D;

    double gamma;
    double gamma_1;
    double cfl;
    double small_pressure;
    double max_dt_increase;
    double max_time;
    double small_dt;

    /*
    void riemann_solver();
    void reconstruction();
    void time_incriment();
    */
    
private:
    /*
    void assign_riemann_solver(Parameters parameter);
    void assign_reconstruction(Parameters parameter);
    void assign_time_stepping(Parameters parameter);
    */
};
