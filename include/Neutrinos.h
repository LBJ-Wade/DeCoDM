#ifndef NEUTRINOS_H
#define NEUTRINOS_H

#include <iostream>
#include <string>
#include <math.h>

#ifndef PI
  #define PI 3.14159265358
#endif

double NeutrinoRate(double E, void* params);
double E_min_v(double E, double m_n);
double NeutrinoFlux(double E_v);
double NeutrinoIntegrand(double E, void* params);
double LoadFluxTable();
double ClearFluxTable();





#endif 