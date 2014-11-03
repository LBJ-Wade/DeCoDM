#include "gsl/gsl_sf_legendre.h"
#include "shared.h"


double ChebyshevP(int order, double x)
{
  switch(order)
    {
    case 0:
      return 1.0;
    case 1:
      return x;
    case 2:
      return 2.0*x*x - 1.0;
    case 3:
      return 4.0*pow(x,3.0) -3.0*x;
    case 4:
      return 8.0*pow(x,4.0) - 8.0*x*x + 1.0;
    case 5:
      return 16.0*pow(x,5.0) - 20.0*pow(x, 3.0) + 5.0*x;
    case 6:
      return 32.0*pow(x,6.0) - 48.0*pow(x,4.0) + 18.0*x*x -1.0;
    case 7:
      return 64.0*pow(x,7.0) - 112.0*pow(x,5.0) + 56.0*pow(x,3.0) - 7.0*x;
    case 8:
      return 128.0*pow(x,8.0) - 256.0*pow(x, 6.0) + 160.0*pow(x, 4.0) - 32.0*x*x + 1.0;
    case 9:
      return 256.0*pow(x,9.0) - 576.0*pow(x,7.0) + 432.0*pow(x,5.0) - 120*pow(x,3.0) + 9.0*x;
    case 10:
      return 512.0*pow(x,10.0) - 1280.0*pow(x,8.0) + 1120.0*pow(x, 6.0) - 400.0*pow(x, 4.0) + 50.0*x*x - 1.0;
    }

  printf("Error in polyf.c - Chebyshev polynomials not yet supported for l > 10...exiting\n");
  exit(0);
}



double polyf(double v, void* params)
{

  double alpha = (v/1000.0);

 // double alpha = v/v_max;

  double logf = 0;

  int i;

  for (i = 0; i < N_terms; i++)
  {
    //std::cout << ((double*)params)[i] << std::endl;
   //logf -= pow(alpha,i)*((double*)params)[i];
   //logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*((double*)params)[i];
   logf -= ChebyshevP(i, 2*alpha-1)*((double*)params)[i];
  }
  return v*v*exp(logf);
}

double polyf1(double v, void* params)
{

  double alpha = (v/1000.0);

  // double alpha = v/v_max;

  double logf = 0;

  int i;

  for (i = 0; i < N_terms; i++)
    {
      //std::cout << ((double*)params)[i] << std::endl;
      //logf -= pow(alpha,i)*((double*)params)[i];
      //logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*((double*)params)[i];
      logf -= ChebyshevP(i, 2*alpha-1)*((double*)params)[i];
    }
  return exp(logf);
}

