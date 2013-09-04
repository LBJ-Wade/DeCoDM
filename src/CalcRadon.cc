#include "RadonTransform.h"



int main ()
{
  //Define dimensions of f
  int N = 64;

  int Nx = N;
  int Ny = N;
  int Nz = N;

  int Ntheta = Ny;
  int Nphi = Nz;

  int N_tot = Nx*Ny*Nz;

  //Define input and output grid
  double* f_in;
  f_in = new double [N_tot];

  double* f_out;
  f_out = new double [N_tot];

  //Set values in the input grid
  assDistribution(f_in, Nx, Ny, Nz);

  //Radon transform
  imageRadonTransform3D(f_out,f_in, Nx, Ny, Nz, Ntheta, Nphi);

  //Add up along phi...


  //Print to file
  printMatrix(f_in, "RadonInput.txt", Nx, Ny, Nz);
  printMatrix(f_out, "RadonOutput.txt", Nx, Ny, Nz);

  //Free up memory
  delete [] f_in;
  delete [] f_out;

}