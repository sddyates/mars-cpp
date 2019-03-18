
#pragma once

#include "data_column.hpp"

void flux_tensor(dataColumn&, const Algorithm&, int, int, int);
void face_flux(dataColumn&, const Grid&, const Algorithm&, int, int, int);
double**** RHSOperator(double****, const Grid&, const Algorithm&);
