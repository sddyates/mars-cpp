
#pragma once

#include "grid.hpp"

class dataColumn{
private:
    double* MakeArray1D(int nx);
    double** MakeArray2D(int nvar, int nx);
public:
    dataColumn(int gz, int nvar, int nx);
    int end, beg;
    int nx_tot, nvar, gz;
    double** U;
    double** V;
    double** flux;
    double** FL;
    double** FR;
    double** UL;
    double** UR;
    double** VL;
    double** VR;
    double* SL;
    double* SR;
    double* pres;
    double* cmax;
};
