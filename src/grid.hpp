
#pragma once

#include <ostream>
#include "parameters.hpp"

class Grid {

    public:

        Grid(const Parameters&);

        // Public functions

        // Public members
        int gz;
        int nvar;

        double* x1;
        double* x2;
        double* x3;

        double* x1_verts;
        double* x2_verts;
        double* x3_verts;

        double x1min;
        double x1max;
        double x2min;
        double x2max;
        double x3min;
        double x3max;

        int nx1;
        int nx2;
        int nx3;

        int nx1_tot;
        int nx2_tot;
        int nx3_tot;

        double dx1;
        double dx2;
        double dx3;

        double dxi[3];
        double dxi_min;

        double dv;
        double da;

        int ibeg;
        int iend;
        int jbeg;
        int jend;
        int kbeg;
        int kend;

        int lower_bc_ibeg;
        int lower_bc_iend;
        int lower_bc_jbeg;
        int lower_bc_jend;
        int lower_bc_kbeg;
        int lower_bc_kend;

        int upper_bc_ibeg;
        int upper_bc_iend;
        int upper_bc_jbeg;
        int upper_bc_jend;
        int upper_bc_kbeg;
        int upper_bc_kend;

        int imax;
        int jmax;
        int kmax;

        void apply_boundaries(double);
        void boundary(double ****, Parameters&);

    private:

        int i;

        void cell_centers(double *, double, double, double, int);
        void cell_edges(double *, double *, double, int);

        double state_vector();
        double construct_flux_column();

        void lower_x1_bc();
        void lower_x2_bc();
        void lower_x3_bc();
        void upper_x1_bc();
        void upper_x2_bc();
        void upper_x3_bc();

};
