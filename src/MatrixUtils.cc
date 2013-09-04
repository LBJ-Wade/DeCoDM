#include "MatrixUtils.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//--------Function Definitions-----------

int getIndex(int i, int j, int k, int Nx, int Ny)
{
 return (i + j*Nx + k*Nx*Ny);
}

int getIndexRowMajor(int i, int j, int k, int Ny, int Nz)
{
 return (k + j*Nz + i*Ny*Nz);
}


void normaliseMatrix(double* f, int Nx, int Ny, int Nz)
{
    double Norm = 0;

    for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
      {
	for (int k = 0; k < Nz; k++)
	{
	  Norm += f[getIndex(i,j,k,Nx,Ny)];
	}
      }
    }

 for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
      {
	for (int k = 0; k < Nz; k++)
	{
	  f[getIndex(i,j,k,Nx,Ny)] /= Norm;
	}
      }
    }

}


void printMatrix(double* f, std::string filename, int Nx, int Ny, int Nz) //Extend to 3-D
{
    std::ofstream outputfile;
    outputfile.open (filename.c_str());
    for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
      {
	for (int k = 0; k < Nz; k++)
	{
	  outputfile << f[getIndex(i,j,k,Nx,Ny)] << "\t";
	}
	outputfile << "\n";
      }
      outputfile << "\n";
    }
    outputfile.close();

}

