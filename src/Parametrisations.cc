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

static int Nvbins;
static double dv;
static std::vector<std::vector<std::vector<double> > > velInt_disc;


struct angparams
{double vmin; int j; Astrophysics* astro;
};

double f0_binned(double v,std::vector<double> g,double* v_edges, int N_bins)
{

  double f = 0;
  for (int i = 0; i < N_bins; i++)
  {
	//std::cout << " " << v_edges[i] <<  "\t" << g[i] << std::endl;
    if ((v >= v_edges[i])&&(v < v_edges[i+1]))
    {
      f += 0.5*g[i]*(pow(v_edges[i+1],2) - pow(v,2));
    }
    else if (v < v_edges[i+1])
    {
      f += 0.5*g[i]*(pow(v_edges[i+1],2) - pow(v_edges[i],2));
    }
  }
  //std::cout << " f = " << f << std::endl;
 return 2*PI*f;
}

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

double signum(double x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;	
}

double J(double beta, double x, double y)
{
	double B = beta;
	double t = (1.0-x*x)*(1.0-B*B) - pow((y-B*x),2.0);
	double A1,A2,A3,A4;
	//if ( y < 1e-10) std::cout << y << std::endl;
	if (t < 1e-15)
	{
		A1 = x*(PI/2.0)*signum((y-B*x));
		A2 = y*(PI/2.0)*signum((x-y*B));
		A3 = (PI/4.0)*signum((1.0-y*y - B*B - x + y*B*(1.0+x))/(y-B));
		A4 = -(PI/4.0)*signum((1.0-y*y - B*B + x - y*B*(1.0-x))/(y+B));
		
	}
	else
	{
		double sqrtt = sqrt(t);
	
		A1 = x*asin((y-B*x)/(sqrt((1.0-x*x)*(1.0-B*B))));
		A2 = y*atan((x-B*y)/sqrtt);
		A3 = 0.5*atan((1.0-y*y - B*B - x + y*B*(1.0+x))/((y-B)*sqrtt));
		A4 = -0.5*atan((1.0-y*y - B*B + x - y*B*(1.0-x))/((y+B)*sqrtt));
		

	}
	/*
	if ((y - B) < 1e-10)
	{
		A3 = (PI/4.0)*signum((1.0-y*y - B*B - x + y*B*(1.0+x))/(y-B));
	}
	if ((y+ B) < 1e-10)
	{
		A4 = -(PI/4.0)*signum((1.0-y*y - B*B + x - y*B*(1.0-x))/(y+B));
	}*/
	
	return A1+A2+A3+A4;

}

//Need to define N_ang somewhere...
//Check factors of 2 and pi
double fintegrand(double v, void * params)
{	
	//Need to change notation to match the other one!!!
	struct angparams * myparams = (struct angparams *)params;
	
	int j = myparams->j + 1;
	Astrophysics* astro = myparams->astro;
	double vmin = myparams->vmin;

	//Need to read in and parse params to load in astro stuff...	
	double result = 0;
	
	//Calculate beta and gamma
	double B = vmin/v;
	double gamma = acos(B);
	int l = ceil(0.99999*2.0*N_ang*gamma/PI);
	
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


double fintegrand_disc(double v, void* params)
{	
	
	//struct angparams * myparams = (struct angparams *)params;
	
	int j = j_bin + 1;
	//int k = k_bin + 1;
	double vmin = *((double*)params);
	
	
	if (v < vmin) return 0;
	//
	if (vmin >= 1000) return 0;
    //vmin += 1e-12; //Prevent some errors with vmin = 0
	

	//Need to read in and parse params to load in astro stuff...	
	double result = 0;
	
	//Calculate beta and gamma
	double B = vmin*1.0/v;
	/*
	double res = 0;
	
	
	double A1 = 0;
	if (B < 1e-15) 
	{
		A1 = PI/2.0;
	}
	else
	{
		A1 = atan(sqrt(1- B*B)/B);
	}
	
	if (j == 1)
	{
		res = PI*fdisc(v,0) + A1*(fdisc(v, 1) - fdisc(v,0));
	}
	if (j == 2)
	{
		res = PI*fdisc(v,1) + A1*(fdisc(v, 0) - fdisc(v,1));
	}
	return 4*PI*v*res;
	*/
	
	double gamma = acos(B);
	int l = ceil(0.999999*2*N_ang*gamma/PI);
	
	//Assume l is odd
	int odd = 1;
	int n = (l+1)/2;

	//std::cout << l << std::endl;
	//std::cout <<  v << " " << gamma << " " << l << " " << 2.0*N_ang*gamma/PI << " " << n <<   std::endl;
	//Check if l is even
	if (l % 2 == 0)
	{
		//std::cout << l << std::endl;
		odd = 0;
		n = l/2;
	}
	
	//std::cout << l << " " << n << std::endl;
	if ((j == 1)&&(k_bin == -1))
	{
		std::cout << " Parameters..." << std::endl;
		std::cout << " B = " << B << std::endl;
		std::cout << " gamma = " << gamma << std::endl;
		std::cout << " j = " << j << std::endl;
		std::cout << " v = " << v << std::endl;
		std::cout << " vmin = " << vmin << std::endl;
		std::cout << " l = " << l << std::endl;
		std::cout << " n = " << n << std::endl;
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
	
	//std::cout << "k_+ = " << k_plus[0] << " " << k_plus[1] << " " << k_plus[2] << std::endl; 
	//std::cout << "k_- = " << k_minus[0] << " " << k_minus[1] << " " << k_minus[2] << std::endl; 
	
	//std::cout << " bin edges: " << bin_edges[0] << " " << bin_edges[1] << " " << bin_edges[2] << " " << bin_edges[3] << std::endl; 
	
	for (int i = 0; i < 3; i++)
	{
		result += fdisc(v, k_plus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		
        if (k_plus[i] != k_minus[i]){
			//std::cout << "How many times is this evaluated!" << std::endl;
			///MARK-123 - check minus signs here...
            result -= fdisc(v,k_plus[i]-1)*(J(B, bin_edges[i+1], cos(PI*k_plus[i]*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*k_plus[i]*1.0/N_ang)));
		}
		
	    result += fdisc(v, k_minus[i]-1)*(0.5*PI*(bin_edges[i+1] - bin_edges[i]));
		if (k_plus[i] != k_minus[i])
		{
			///MARK-123 - check minus signs here...
			result += fdisc(v, k_minus[i]-1)*(+J(B, bin_edges[i+1], cos(PI*(k_minus[i]-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k_minus[i]-1)*1.0/N_ang)));
		}
		//Check that this is equivalent to the python...
		for (int k = (k_plus[i]+1); k < k_minus[i]; k++)
		{
			result += fdisc(v,k-1)*(+ J(B, bin_edges[i+1], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i], cos(PI*(k-1)*1.0/N_ang)) - J(B, bin_edges[i+1], cos(PI*(k)*1.0/N_ang)) + J(B, bin_edges[i], cos(PI*(k)*1.0/N_ang)));
		}
	}
 	
	return 4.0*PI*v*result;
}

double velInt_DRT(double vmin,  Astrophysics* astro)
{
	//Need to very carefully implement the integral numerically...
	
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
	
	//MARK-123
	//Should this be 'int n = 0'?
	for (int n = 1; n < N_ang+1; n++)
	{
		//Sort all this out...
		//Print values...
		
		v1 = std::min(vmin/cos((n-1)*PI*0.5/N_ang),1000.0);
		v2 = std::min(vmin/cos((n)*PI*0.5/N_ang), 1000.0);
    	status = gsl_integration_qag(&F,v1,v2, 0.1, 1e-3, 5000,4,
				workspace, &tempres, &error);
				result += tempres;
	}
	return result;		  
}

double velInt_DRT_disc(double vmin,  Astrophysics* astro)
{
	
	
	static int count = 0;
    if (count == 0)
	{
		int tempj = j_bin;
		//std::cout << " Edit: MARK-123" << std::endl;
		std::cout << " Note: I've changed the normalisation to be 3 separate evals..." << std::endl;
		std::cout << " Check that Nvbins matches that in CalcLikelihood.cc" << std::endl;
		Nvbins = 100;
		dv = 1000.0/Nvbins;
		calcApproxMatrix();
		count++;
		j_bin = tempj;
	}

	double fudge = 1.0;
	int vqi = floor(vmin/dv);
	//std::cout << vmin << "\t" << vqi << std::endl;
	if (vqi >= (Nvbins - 1)) return 0;
	//std::cout << vmin << "\t" << vqi << std::endl;
	double result = 0;
	for (int k = 0; k < N_ang; k++)
	{
		int jkind = N_ang*j_bin + k;
		//std::cout << jkind << std::endl;
		//for (int vi = vqi; vi < Nvbins; vi++)
		for (int vi = 0; vi < Nvbins; vi++)
		{
			if (k == 1)
			{
				fudge = 1.0;
			}
			//Could replace polyf_ang with an integrated version...
			//I should probably do that...
			//Linear interpolation...?
			//result += (velInt_disc[jkind][vi][vqi] + (velInt_disc[jkind][vi][vqi+1]- velInt_disc[jkind][vi][vqi])*(vmin/dv - vqi))*polyf_ang(dv*vi,astro,k);
			//Check this linear interpolation...
			//if (j_bin == 1) fudge = 0;
			result += fudge*(velInt_disc[jkind][vi][vqi] + (velInt_disc[jkind][vi][vqi+1] - velInt_disc[jkind][vi][vqi])*(vmin/dv - vqi))*astro->bin_params[k][vi];
			//result += velInt_disc[jkind][vi][vqi]*astro->bin_params[k][vi];
		}		
	}
	return result;
}


//N_ang x N_ang matrix, each element is a vector of 1000 elements, 
//corresponding to bins in v...

double fdisc(double v, int k)
{
	//Should this be k_bin or k_bin+1

	if (k != k_bin) return 0;
	//Check that this works correctly with rounding...(ceil or floor?)
	//std::cout << floor(v/dv) << " " << ind << std::endl;
	if (floor(v/dv) != ind) return 0; //1000 bins in v from 0 to 1000...
	return 1.0;
}


double calcApproxMatrix()
{
	//Code in a specific parameter for number of bins...
	
	//Initialise the integrator
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
    F.function = &fintegrand_disc;

    result = 0;
	
	std::vector<double> velInt1;
	std::vector<std::vector<double> > velInt2;
	
	std::cout << " Calculating approximate angular velocity integral..." << std::endl;
	velInt_disc.clear();
	for (int j = 0; j < N_ang; j++)
	{
		j_bin = j;
		for (int k = 0; k < N_ang; k++)
		{
			std::cout << "\t j = " << j+1 << "; k = " << k+1 << std::endl;
			k_bin = k;
			velInt2.clear();
			for (int vi = 0; vi < Nvbins; vi++)
			{	
				ind = vi;
				
				velInt1.clear();
				for (int vqi = 0; vqi < Nvbins; vqi++)
				{
					double vmin = dv*vqi;
					//std::cout << vmin << std::endl;
					F.params = &vmin;
					//vmin = 100;
					//std::cout << fintegrand_disc(200, &vmin) << std::endl;
					result = 0;
					//Only need to integrate over bin that's non-zero
			    	//status = gsl_integration_qag(&F,vmin,1000.0, 0.1, 1e-6, 5000,4,
					
					status = gsl_integration_qag(&F,vi*dv,(vi+1.0)*dv, 0, 1e-6, 5000,4,
							workspace, &result, &error);
					
							velInt1.push_back(result);
					
					//velInt1.push_back(fintegrand_disc(dv*vi, &vmin));
							
							//std::cout << vi*dv << " " << (vi+1.0)*dv << std::endl;
							//std::cout << result << std::endl;
				}
				velInt2.push_back(velInt1);
			}
			velInt_disc.push_back(velInt2);
		}
	}	
	std::cout << " Finished calculating." << std::endl;
}


#endif