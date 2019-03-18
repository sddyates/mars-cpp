
#pragma once

inline double**** MakeArray4D(int nvar, int nx3, int nx2, int nx1) {

    /*
    int var, k, j;

    Array4D = new double[nvar];

    for (int var = 0; var < nvar; var++) {
        Array4D[var] = new double**[nx3];

          for (int k = 0; k < nx3; ++k) {
              Array4D[var][k] = new double*[nx2];

              for (int j = 0; j < nx3; ++j) {
                  Array4D[var][k][j] = new double[nx1];
              }
          }
    }

    return Array4D;
    */

    int var, k, j;
    double**** Array4D = new double***[nvar];
    for (var = 0; var < nvar; var++) {
        Array4D[var] = new double**[nx3];

          for (k = 0; k < nx3; ++k) {
              Array4D[var][k] = new double*[nx2];

              for (j = 0; j < nx2; ++j) {
                  Array4D[var][k][j] = new double[nx1];
              }
          }
    }

    return Array4D;
}

inline double** MakeArray2D(int nvar, int nx) {

    int i;
    double** Array2D = new double*[nvar];
    for (i = 0; i < nvar; ++i) {
        Array2D[i] = new double[nx];
    }

    return Array2D;
}

inline double* MakeArray1D(int nvar) {

    double* Array1D = new double[nvar];

    return Array1D;
}

// Generate state vector to hold conservative
// and primative variables.
//double V[grid.nvar][grid.nx3][grid.nx2][grid.nx1];
//vector<vector<vector<vector<double> > > > V;

// Set up sizes. (LENGTH x HEIGHT x WIDTH x DEPTH)
/*V.resize(grid.nvar);
for (int var = 0; var < grid.nvar; ++var) {
    V[var].resize(grid.nx3);

    for (int k = 0; k < grid.nx3; ++k) {
        V[var][k].resize(grid.nx2);

        for (int j = 0; j < grid.nx2; ++j) {
            V[var][k][j].resize(grid.nx1);
        }
    }
}*/
