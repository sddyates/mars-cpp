
#include <time.hpp>

void Time::Time () {
    t = 0.0:
    time->dt = parameter.initial_dt;
    i = 0;
    num = 1;
    percent = 100.0/problem.parameter.max_time;
}

void Time::iterate () {

}

void Time::compare () {
    t = 0.0:
    time->dt = parameter.initial_dt;
    i = 0;
    num = 1;
    percent = 100.0/problem.parameter.max_time;


    dt = min(dt_new, problem.parameter.max_dt_increase*dt);

    if ((t + dt) > problem.parameter.max_time) {
        dt = problem.parameter.max_time - t;
    }

    if ((t + dt) > num*problem.parameter.plot_frequency) {
        dump(V, grid, problem.parameter, num);
        num += 1
    }

    t += dt;
    i++;
}

void time_step(double V, Grid &g, const Algorithm &a, TimeStructure &time) {

    double cs, cs_max = 0.0;
    double velocity, velocity_max = 0.0;
    double speed_max;
    int i, j, k;

    for (k=g.kbeg; k=g.kend; k++) {
        for (j=g.jbeg; j=g.jend; j++) {
            for (i=g.ibeg; i=g.iend; i++) {

                cs = std::sqrt(a.gamma*V[prs][k][j][i]/V[rho][k][j][i]);

                if (cs > cs_max) {
                    cs_max = cs;
                }

                velocity = std::max(fabs(V[vx1][k][j][i]), fabs(V[vx2][k][j][i]))
                velocity = std::max(fabs(velocity), fabs(V[vx3][k][j][i]))

                if (velocity > velocity_max) {
                    velocity_max = velocity;
                }

            }
        }
    }

    speed_max = velocity_max + cs_max;
    dt = a.cfl*g.min_dxi/speed_max;
    mach_number = velocity_max/cs_max;

    if (dt < small_dt || isnan(dt)) {
        std::cout << "Error, dt is nan to small, dt = " << dt << std::endl;
    }

    time->dt = dt;
    time->velocity_max = velocity_max;
    time->mach_number = mach_number;

    return;
}
