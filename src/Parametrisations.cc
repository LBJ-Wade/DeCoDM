#ifndef PARAMETRISATIONS_CC
#define PARAMETRISATIONS_CC

#include "Parametrisations.h"
#include <math.h>
#include <iostream>

#ifndef PI
  #define PI 3.14159265358
#endif

double f0_binned(double v,double* g,double* v_edges, int N_bins)
{

  double f = 0;
  for (int i = 0; i < N_bins; i++)
  {
    if ((v >= v_edges[i])&&(v < v_edges[i+1]))
    {
      f += 0.5*g[i]*(pow(v_edges[i+1],2) - pow(v,2));
    }
    else if (v < v_edges[i+1])
    {
      f += 0.5*g[i]*(pow(v_edges[i+1],2) - pow(v_edges[i],2));
    }
  }

 return 2*PI*f;
}

double eta_binned(double v,double* g,double* v_edges, int N_bins)
{
  double f = 0;
  for (int i = 0; i < N_bins; i++)
  {
    if ((v >= v_edges[i])&&(v < v_edges[i+1]))
    {
      f = g[i];
    }
  }

 return f;
}

double f1_binned(double v,double* g,double* v_edges, int N_bins)
{
  double f = 0;
  for (int i = 0; i < N_bins; i++)
  {
    if ((v >= v_edges[i])&&(v < v_edges[i+1]))
    {
      f += g[i]*(v_edges[i+1] - v);
    }
    else if (v < v_edges[i+1])
    {
      f += g[i]*(v_edges[i+1] - v_edges[i]);
    }
  }

 return 2*PI*v*f;

}

double f2_binned(double v,double* g,double* v_edges, int N_bins)
{

  double f = 0;
  for (int i = 0; i < N_bins; i++)
  {
    if ((v >= v_edges[i])&&(v < v_edges[i+1]))
    {
      f += g[i]*log(v_edges[i+1]/v);
    }
    else if (v < v_edges[i+1])
    {
      f += g[i]*log(v_edges[i+1]/v_edges[i]);
    }
  }

 return 3*PI*v*v*f + 0.5*f0_binned(v,g,v_edges,N_bins);
}


#endif