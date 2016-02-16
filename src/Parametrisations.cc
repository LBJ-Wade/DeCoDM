#ifndef PARAMETRISATIONS_CC
#define PARAMETRISATIONS_CC

#include "DMUtils.h"
#include "Parametrisations.h"
#include "Distributions.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include <math.h>
#include <iostream>

#ifndef PI
  #define PI 3.14159265358
#endif

int ind;
int k_bin;

struct angparams
{double vmin; int j; Astrophysics* astro;
};

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

double J(double beta, double x, double y)
{
	double B = beta;
	double sqrtt = sqrt((1.0-x*x)*(1.0-B*B) - pow((y-B*x),2.0));
	
	double A1 = x*asin((y-B*x)/(sqrt((1.0-x*x)*(1.0-B*B))));
	double A2 = y*atan((x-B*y)/sqrtt);
	double A3 = 0.5*atan((1.0-y*y - B*B - x + y*B*(1.0+x))/((y-B)*sqrtt));
	double A4 = -0.5*atan((1.0-y*y - B*B + x - y*B*(1.0-x))/((y+B)*sqrtt));
	return A1+A2+A3+A4;
}

//Need to define N_ang somewhere...
//Check factors of 2 and pi
double fintegrand(double v, void * params)
{	
	
	struct angparams * myparams = (struct angparams *)params;
	
	int j = myparams->j + 1;
	Astrophysics* astro = myparams->astro;
	double vmin = myparams->vmin;

	//Need to read in and parse params to load in astro stuff...	
	double result = 0;
	
	//Calculate beta and gamma
	double B = vmin/v;
	double gamma = acos(B);
	int l = ceil(2.0*N_ang*gamma/PI);
	
	//Assume l is odd
	int odd = 1;
	int n = (l+1)/2;

	//Check if l is even
	if (l % 2 == 0)
	{
		odd = 0;
		n = l/2;
	}
	
	//Initialise the k-values
	int k_plus[3];
	int k_minus[3];
	
	//Note that the order of these is reversed w.r.t.
	//the python version
	if (j <= (n-1))
	{
        k_plus[0] = n-j;
        k_plus[1] = n-j+1 - odd;				
        k_plus[2] = n-j+1;
	}
	else if (j == n)
	{
        k_plus[0] = 1;
        k_plus[1] = 1;				
        k_plus[2] = 1;
	}
	else
	{
        k_plus[0] = j-n+1;
        k_plus[1] = j-n + odd;			
        k_plus[2] = j-n;
	}
	
    if (j <= (N_ang-n))
	{
        k_minus[0] = n+j;		
        k_minus[1] = n+j - odd;
        k_minus[2] = n+j - 1;
	}
    else if (j == (N_ang+1-n))
	{
        k_minus[0] = N_ang;
	    k_minus[1] = N_ang;		
        k_minus[2] = N_ang;
	}
    else
	{
		k_minus[0] = 2*N_ang-n-j+1;
	    k_minus[1] = 2*N_ang-n-j+1 + odd;
        k_minus[2] = 2*N_ang-n-j+2;
	}
	
	//Initialise the bin edges
	double bin_edges[4];
	bin_edges[0] = cos(j*PI/N_ang);
	bin_edges[3] = cos((j-1.0)*PI/N_ang);
	
    if (odd)
	{
        bin_edges[1] = cos((j + n - 1)*PI/N_ang - gamma);
        bin_edges[2] = cos((j - n)*PI/N_ang + gamma);
	}    
    else
	{
        bin_edges[1] = cos((j - n)*PI/N_ang + gamma);
        bin_edges[2] = cos((j + n - 1)*PI/N_ang - gamma);
	}
	
	for (int i = 0; i < 3; i++)
	{
		result += (2.0/v)*polyf_ang(v, astro, k_plus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		
        if (k_plus[i] != k_minus[i]){
            result -= (2.0/v)*polyf_ang(v, astro, k_plus[i]-1)*(J(B, bin_edges[i+1], cos(PI*k_plus[i]*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*k_plus[i]*1.0/N_ang)));
		}
		
	    result += (2.0/v)* polyf_ang(v, astro, k_minus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		if (k_plus[i] != k_minus[i])
		{
			result += (2.0/v)*polyf_ang(v, astro, k_minus[i]-1)*(+J(B, bin_edges[i+1], cos(PI*(k_minus[i]-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k_minus[i]-1)*1.0/N_ang)));
		}
		
		//Check that this is equivalent to the python...
		for (int k = (k_plus[i]+1); k < k_minus[i]; k++)
		{
			result += (2.0/v)*polyf_ang(v, astro,k-1)*(+ J(B, bin_edges[i+1], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i+1], cos(PI*(k)*1.0/N_ang)) + J(B, bin_edges[i], cos(PI*(k)*1.0/N_ang)));
		}
	}
 	
	return v*v*result*2.0*PI;
}

double fintegrand_disc(double v, void * params)
{	
	
	struct angparams * myparams = (struct angparams *)params;
	
	int j = j_bin + 1;
	//int k = k_bin + 1;
	double vmin = *((double*)params);

	//Need to read in and parse params to load in astro stuff...	
	double result = 0;
	
	//Calculate beta and gamma
	double B = vmin/v;
	double gamma = acos(B);
	int l = ceil(2.0*N_ang*gamma/PI);
	
	//Assume l is odd
	int odd = 1;
	int n = (l+1)/2;

	//Check if l is even
	if (l % 2 == 0)
	{
		odd = 0;
		n = l/2;
	}
	
	//Initialise the k-values
	int k_plus[3];
	int k_minus[3];
	
	//Note that the order of these is reversed w.r.t.
	//the python version
	if (j <= (n-1))
	{
        k_plus[0] = n-j;
        k_plus[1] = n-j+1 - odd;				
        k_plus[2] = n-j+1;
	}
	else if (j == n)
	{
        k_plus[0] = 1;
        k_plus[1] = 1;				
        k_plus[2] = 1;
	}
	else
	{
        k_plus[0] = j-n+1;
        k_plus[1] = j-n + odd;			
        k_plus[2] = j-n;
	}
	
    if (j <= (N_ang-n))
	{
        k_minus[0] = n+j;		
        k_minus[1] = n+j - odd;
        k_minus[2] = n+j - 1;
	}
    else if (j == (N_ang+1-n))
	{
        k_minus[0] = N_ang;
	    k_minus[1] = N_ang;		
        k_minus[2] = N_ang;
	}
    else
	{
		k_minus[0] = 2*N_ang-n-j+1;
	    k_minus[1] = 2*N_ang-n-j+1 + odd;
        k_minus[2] = 2*N_ang-n-j+2;
	}
	
	//Initialise the bin edges
	double bin_edges[4];
	bin_edges[0] = cos(j*PI/N_ang);
	bin_edges[3] = cos((j-1.0)*PI/N_ang);
	
    if (odd)
	{
        bin_edges[1] = cos((j + n - 1)*PI/N_ang - gamma);
        bin_edges[2] = cos((j - n)*PI/N_ang + gamma);
	}    
    else
	{
        bin_edges[1] = cos((j - n)*PI/N_ang + gamma);
        bin_edges[2] = cos((j + n - 1)*PI/N_ang - gamma);
	}
	
	for (int i = 0; i < 3; i++)
	{
		result += (2.0/v)*fdisc(v, k_plus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		
        if (k_plus[i] != k_minus[i]){
            result -= (2.0/v)*fdisc(v,k_plus[i]-1)*(J(B, bin_edges[i+1], cos(PI*k_plus[i]*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*k_plus[i]*1.0/N_ang)));
		}
		
	    result += (2.0/v)*fdisc(v, k_minus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		if (k_plus[i] != k_minus[i])
		{
			result += (2.0/v)*fdisc(v, k_minus[i]-1)*(+J(B, bin_edges[i+1], cos(PI*(k_minus[i]-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k_minus[i]-1)*1.0/N_ang)));
		}
		
		//Check that this is equivalent to the python...
		for (int k = (k_plus[i]+1); k < k_minus[i]; k++)
		{
			result += (2.0/v)*fdisc(v,k-1)*(+ J(B, bin_edges[i+1], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i+1], cos(PI*(k)*1.0/N_ang)) + J(B, bin_edges[i], cos(PI*(k)*1.0/N_ang)));
		}
	}
 	
	return v*v*result*2.0*PI;
}

double velInt_DRT(double vmin,  Astrophysics* astro)
{
	struct angparams params = {vmin, j_bin, astro};
	
    //Declare gsl workspace (1000 subintervals)
    gsl_integration_workspace * workspace
           = gsl_integration_workspace_alloc (5000);
    //switch off default error handler, store old error handler in
    //old_handler
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    //Declare gsl function to be integrated
    gsl_function F;

    double result, tempres, error;
    int status;
	
	//auto fj = [=] (double v, void* params){ return fintegrand(v, params, vmin, j);};

	//Need to pass something through paramvals...

	F.params = &params;
    F.function = &fintegrand;

    result = 0;
	double v1,v2;
	for (int n = 0; n < N_ang+1; n++)
	{
		v1 = std::min(vmin/cos((n-1)*PI*0.5/N_ang),1000.0);
		v2 = std::min(vmin/cos((n)*PI*0.5/N_ang), 1000.0);
    	status = gsl_integration_qag(&F,v1,v2, 0.1, 1e-3, 5000,4,
                               workspace, &tempres, &error);
		result += tempres;
	}
	return result;		  
}


//N_ang x N_ang matrix, each element is a vector of 1000 elements, 
//corresponding to bins in v...

double fdisc(v, k)
{
	if (k != k_bin) return 0;
	if (floor(v) != ind) return 0; //1000 bins in v from 0 to 1000...
	return 1;
}


double calcApproxMatrix()
{
	std::vector<std::vector<double> > velInt_disc;
	
	
	
}



#endif