#ifndef RADONTRANSFORM_H
#define RADONTRANSFORM_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>
#include <string.h>
#include "MatrixUtils.h"

#ifndef PI
  #define PI 3.14159265358
#endif


//Working reasonably well - only shows a weak dependence on the number of pixels...

//At the moment everything is in degrees!!!
//i,j are zero based indices
//Currently, q is an 'integer' index...change so that I specify at different values...

//Pixel inversion...

//Event generation

//Event rate calculation...

//Maybe do some interpolation so that I can specify Nq...

//Maybe redo to allow in-place rotations and transforms...

//Use vectors...


//--------Function Prototypes----------

void rotateImage(double* Rout, double* Rin, double theta, int Nx, int Ny, int Nz, int axis);

double rotatePointX(int i, int j, double theta);
double rotatePointY(int i, int j, double theta);
double rotatePointZ(int i, int j, double theta);

double pointRadonTransform(double* fin, int Nx, int Ny, int Nz, double q, double theta, double phi);

void imageRadonTransform(double* fout, double* fin, int Nx, int Ny, int Nz);

void imageRadonTransform2D(double* fout, double* fin, int Nx, int Ny, int Ntheta);
void imageRadonTransform3D(double* fout, double* fin, int Nx, int Ny, int Nz, int Ntheta, int Nphi);

void inverseRadonTransform2D(double* fout, double* fin, int Nx, int Ny, int Ntheta);
void inverseRadonTransform3D(double* fout, double* fin, int Nx, int Ny, int Nz, int Ntheta, int Nphi);

void applyFilter2D(double* fout, double* fin, int Nx, int Ny);
void applyFilter3D(double* fout, double* fin, int Nx, int Ny, int Nz);

void Laplacian (double* fout, double* fin, int Nx, int Ny, int Nz);

double mass(double* f, int Nq, int Ntheta, int Nphi, int theta, int phi);
double CoM(double* f, int Nq, int Ntheta, int Nphi, int theta, int phi);

void maxwellDistribution(double* f, int Nx, int Ny, int Nz, double vlagX, double vlagY, double vlagZ, double sigma);
void assDistribution(double* f, int Nx, int Ny, int Nz);

#endif