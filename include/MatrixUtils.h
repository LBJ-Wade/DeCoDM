#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>

#ifndef PI
  #define PI 3.14159265358
#endif

//--------Function Prototypes------------
int getIndex(int i, int j, int k, int Nx, int Ny);
int getIndexRowMajor(int i, int j, int k, int Ny, int Nz);

void normaliseMatrix(double* f, int Nx, int Ny, int Nz);
void printMatrix(double* f, std::string filename, int Nx, int Ny, int Nz);

#endif