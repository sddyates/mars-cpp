
#pragma once

struct Parameters {

    const char* name;

    const char* dimensions;
    double x1_min;
    double x1_max;
    double x2_min;
    double x2_max;
    double x3_min;
    double x3_max;
    int x1_resolution;
    int x2_resolution;
    int x3_resolution;

    double cfl;
    double initial_dt;
    double max_dt_increase;
    double max_time;

    double plot_frequency;
    bool print_to_file;

    double gamma;
    double unit_density;
    double unit_length;
    double unit_velocity;

    const char* riemann;
    const char* reconstruction;
    const char* limiter;
    const char* time_stepping;
    const char* physics;

    const char* lower_x1_boundary;
    const char* lower_x2_boundary;
    const char* lower_x3_boundary;
    const char* upper_x1_boundary;
    const char* upper_x2_boundary;
    const char* upper_x3_boundary;
    bool internal_boundary;

};

inline void initilise_parameters(Parameters &parameter) {

    parameter.name = "shock tube";

    parameter.dimensions = "1D";
    parameter.x1_min = 0.0;
    parameter.x1_max = 1.0;
    parameter.x2_min = 0.0;
    parameter.x2_max = 1.0;
    parameter.x3_min = 0.0;
    parameter.x3_max = 0.0;
    parameter.x1_resolution = 8;
    parameter.x2_resolution = 8;
    parameter.x3_resolution = 1;

    parameter.cfl = 0.3;
    parameter.initial_dt = 1.0e-5;
    parameter.max_dt_increase = 1.5;
    parameter.max_time = 1.0;

    parameter.plot_frequency = 0.1;
    parameter.print_to_file = false;

    parameter.gamma = 1.666666;
    parameter.unit_density = 1.0;
    parameter.unit_length = 1.0;
    parameter.unit_velocity = 1.0;

    parameter.riemann = "tvdlf";
    parameter.reconstruction = "linear";
    parameter.limiter = "minmod";
    parameter.time_stepping = "RK2";
    parameter.physics = "hd";

    parameter.lower_x1_boundary = "outflow";
    parameter.lower_x2_boundary = "outflow";
    parameter.lower_x3_boundary = "outflow";
    parameter.upper_x1_boundary = "outflow";
    parameter.upper_x2_boundary = "outflow";
    parameter.upper_x3_boundary = "outflow";

    parameter.internal_boundary = false;

    return;
}
