
#include "data_column.hpp"

dataColumn::dataColumn(int gzi, int nvari, int nxi) {

    U = MakeArray2D(nvari, nxi);
    V = MakeArray2D(nvari, nxi);
    flux = MakeArray2D(nvari, nxi);
    FL = MakeArray2D(nvari, nxi);
    FR = MakeArray2D(nvari, nxi);
    UL = MakeArray2D(nvari, nxi);
    UR = MakeArray2D(nvari, nxi);
    VL = MakeArray2D(nvari, nxi);
    VR = MakeArray2D(nvari, nxi);
    SL = MakeArray1D(nxi);
    SR = MakeArray1D(nxi);
    pres = MakeArray1D(nxi);
    cmax = MakeArray1D(nxi);

    nx_tot = nxi;
    nvar = nvari;
    gz = gzi;

    beg = gzi;
    end = nxi - 2*gzi;
}

double *dataColumn::MakeArray1D(int nx) {

    double* Array1D = new double[nx];

    return Array1D;
}

double **dataColumn::MakeArray2D(int nvar, int nx) {

    int i;
    double** Array2D = new double*[nvar];
    for (i = 0; i < nvar; ++i) {
        Array2D[i] = new double[nx];
    }

    return Array2D;
}
