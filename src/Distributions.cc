
#include "Astrophysics_Class.h"
#include "Distributions.h"

#include "gsl/gsl_sf.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"
#include "DMUtils.h"
#include "Parametrisations.h"

#ifndef PI
  #define PI 3.14159265358
#endif

double v_max = 1000.0;

//-------------------------------------------------------------------------
//-----------Velocity Integrals--------------------------------------------
//-------------------------------------------------------------------------

double velInt_maxwell(double v, Astrophysics* astro)
{


  double vel_integral_total = 0;

  for (int i = 0; i < astro->N_dist; i++)
  {
          double v_rms = astro->v_rms[i];
          double v_esc = astro->v_esc;
          double v_lag = astro->v_lag[i];
		  
          double vel_integral = 0;
          double N = 1.0/(gsl_sf_erf(v_esc/(sqrt(2)*v_rms)) - sqrt(2.0/PI)*(v_esc/v_rms)*exp(-0.5*pow(v_esc/v_rms,2)));
	  /*if (v > (v_esc + v_lag))
	  {
	   return 0;
	  }
	  else if (v < (v_esc - v_lag))
	  {
	    vel_integral = (0.5/(v_lag))*(gsl_sf_erf((v + v_lag)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
	    vel_integral -= sqrt(2.0/PI)*(1.0/v_rms)*exp(-0.5*pow(v_esc/v_rms,2));
	  }
	  else
	  {
	    vel_integral = (0.5/(v_lag))*(gsl_sf_erf((v_esc)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
	    vel_integral -= 0.5*sqrt(2.0/PI)*((v_esc+v_lag - v)/(v_lag*v_rms))*exp(-0.5*pow(v_esc/v_rms,2));
	  }*/
		  
		  double aplus = std::min((v + v_lag), v_esc)/(sqrt(2)*v_rms);
		  double aminus = std::min((v - v_lag), v_esc)/(sqrt(2)*v_rms);

		  vel_integral = (0.5/v_lag)*(gsl_sf_erf(aplus) - gsl_sf_erf(aminus));
		  vel_integral -= (1.0/(sqrt(PI)*v_lag))*(aplus - aminus)*exp(-0.5*pow(v_esc/v_rms, 2));

          vel_integral_total += 2*PI*N*vel_integral*astro->fraction[i];
  }
  return vel_integral_total;
}

double velInt_maxwell_modified(double v, Astrophysics* astro)
{


  double vel_integral_total = 0;

  for (int i = 0; i < astro->N_dist; i++)
  {
          double v_rms = astro->v_rms[i];
          double v_esc = astro->v_esc;
          double v_lag = astro->v_lag[i];

          double vel_integral = 0;
          
		  double T_plus = gsl_sf_erf((v + v_lag)/(sqrt(2)*v_rms));
		  double T_minus = gsl_sf_erf((v - v_lag)/(sqrt(2)*v_rms));
		  
		  double K_plus = exp(-0.5*pow((v+v_lag)/v_rms,2));
		  double K_minus = exp(-0.5*pow((v-v_lag)/v_rms,2));

          vel_integral = (sqrt(PI/2.0)/v_lag)*(sqrt(2.0*PI)*(v_rms*v_rms - v*v + v_lag*v_lag)*(T_plus - T_minus) + 2.0*v_rms*(v_lag - v)*K_plus + 2.0*v_rms*(v+v_lag)*K_minus);

          vel_integral_total += vel_integral*astro->fraction[i]/(3e5*3e5);
  }
  return vel_integral_total;
}

//Returns the unnormalised Lisanti et al directionally averaged speed distribution
double Lisanti_f(double v, void* params)
{
  Astrophysics* astro = ((Astrophysics*)params);

  double v_esc = astro->v_esc;
  double v0 = astro-> v0;
  int k = astro->k;

  if (v > (v0 + v_esc)) return 0;

 double tot = 0;

 double A = v_esc*v_esc - v*v - v0*v0;
 double B = 2*v0*v;
 double C = k*v0*v0;

 if (v < (v_esc - v0))
 {
    for (int i = 0; i < k; i++)
    {
      tot += 2*PI*(C/B)*pow(-1,i)*(nCr(k, i)/(1.0*k-i))*(exp((k-i)*(A+B)/C) - exp((k-i)*(A-B)/C));
    }
    tot += pow(-1,k)*4*PI;

 }
 else
 {
   for (int i = 0; i < k; i++)
    {
      tot += 2*PI*(C/B)*pow(-1,i)*(nCr(k, i)/(1.0*k-i))*(exp((k-i)*(A+B)/C) - 1);
    }
    tot += 2*PI*pow(-1,k)*(1+A/B);
 }

 return v*v*tot;
}

double LisantiIntegrand(double v, void* params)
{
  Astrophysics* astro = ((Astrophysics*)params);

  double v_esc = astro->v_esc;
  double v0 = astro-> v0;
  int k = astro->k;


  if (v > (v0 + v_esc)) return 0;

  double tot = 0;

  double A = v_esc*v_esc - v*v - v0*v0;
  double B = 2*v0*v;
  double C = k*v0*v0;

  if (v < (v_esc - v0))
  {
      for (int i = 0; i < k; i++)
      {
	tot += 2*PI*(C/B)*pow(-1,i)*(nCr(k, i)/(1.0*k-i))*(exp((k-i)*(A+B)/C) - exp((k-i)*(A-B)/C));
      }
      tot += pow(-1,k)*4*PI;

  }
  else
  {
    for (int i = 0; i < k; i++)
      {
	tot += 2*PI*(C/B)*pow(-1,i)*(nCr(k, i)/(1.0*k-i))*(exp((k-i)*(A+B)/C) - 1);
	//std::cout << tot << std::endl;
      }
      tot += 2*PI*pow(-1,k)*(1+A/B);
  }
  return v*tot;
}


//Returns the norm of the Lisanti distribution (i.e. the integral over all speeds)
double Lisanti_norm(void* params)
{

  //Declare gsl workspace (5000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

  //Declare gsl function to be integrated
  gsl_function F;
  F.function = &Lisanti_f;

  F.params = params;

  double result, error;

    int status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 5000,6,
                             workspace, &result, &error);

			     //if (result < 0) std::cout << "Negative rate!" << std::endl;

  if (status ==  GSL_EROUND)
  {
  //result = 0;
  std::cout << "GSL rounding error!" << std::endl;
  std::cout << result << std::endl;
  }

  //Free workspace
  gsl_integration_workspace_free (workspace);

  return result;
}


//Returns the UNNORMALISED rate based on a Lisanti speed distribution
double velInt_Lisanti(double v, Astrophysics* astro)
{
   double vel_integral = integrator(v, &LisantiIntegrand, astro);
   return vel_integral;
}

double velInt_isotropicBinned(double v, Astrophysics* astro)
{
    double vel_integral = f0_binned(v,astro->vel_params,astro->bin_edges, astro->N_vp);
    return vel_integral;
}

double velInt_forwardBinned(double v, Astrophysics* astro)
{
  //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, astro->vel_params_forward, astro->N_vp);
  //And backward rate
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, astro->vel_params_backward, astro->N_vp);
  return vel_integral;
}

double velInt_backwardBinned(double v, Astrophysics* astro)
{
  //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, astro->vel_params_backward, astro->N_vp);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, astro->vel_params_forward, astro->N_vp);
  return vel_integral;
}

double velInt_isotropicPoly(double v, Astrophysics* astro)
{
    double vel_integral = multipoleRadon(v, 0, &multipoleIntegrand,astro->vel_params, astro->N_vp);
    return vel_integral;
}

double velInt_forwardPoly(double v, Astrophysics* astro)
{
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, astro->vel_params_forward, astro->N_vp);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, astro->vel_params_backward, astro->N_vp);
  return vel_integral;
}


double velInt_backwardPoly(double v, Astrophysics* astro )
{
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, astro->vel_params_backward, astro->N_vp);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, astro->vel_params_forward, astro->N_vp);
  return vel_integral;
}

double velInt_polytotal(double v, Astrophysics* astro)
{
    if (v > v_max) return 0;

    //Declare gsl workspace (1000 subintervals)
    gsl_integration_workspace * workspace
           = gsl_integration_workspace_alloc (5000);

    //Declare gsl function to be integrated
    gsl_function F;
    F.function = &polyfintegrand_total;

    F.params = astro;

    double result, error;

    int status = gsl_integration_qag(&F,v,v_max, 0, 1e-6, 5000,6,
                               workspace, &result, &error);

    if (status ==  GSL_EROUND)
    {
    //result = 0;
    std::cout << "GSL rounding error!" << std::endl;
    std::cout << result << std::endl;
    }

    //Free workspace
    gsl_integration_workspace_free (workspace);

    //Return result of integration
    return 2.0*PI*result;
}

double multipoleRadon(double v_q, int l, double integrand (double,void*), double* params, int N)
{
  if (v_q > v_max)
  {
   return 0;
  }

  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

  //Declare gsl function to be integrated
  gsl_function F;
  F.function = integrand;

  //Pack paramaters
  void *full_params[4];

  full_params[0] = malloc(sizeof(int));
  full_params[1] = malloc(sizeof(int));
  full_params[2] = malloc(sizeof(double));
  full_params[3] = malloc(N*sizeof(double));

  *((int*)full_params[0]) = l;
  *((int*)full_params[1]) = N;
  *((double*)full_params[2]) = v_q;
  *((double**)full_params[3]) = params;

  F.params = full_params;

  double result, error;
  //std::cout << "Start" << std::endl;
  //std::cout << v_q << std::endl;

  //switch off default error handler, store old error handler in
  //old_handler
   //gsl_error_handler_t * old_handler=gsl_set_error_handler_off();


    int status = gsl_integration_qag(&F,v_q,v_max, 0, 1e-6, 5000,6,
                             workspace, &result, &error);

			     //if (result < 0) std::cout << "Negative rate!" << std::endl;

  if (status ==  GSL_EROUND)
  {
  //result = 0;
  std::cout << "GSL rounding error!" << std::endl;
  std::cout << result << std::endl;
  }

  //Free workspace
  gsl_integration_workspace_free (workspace);

  //////////////////////-------------------------Do I need to free more of these...?

  free(full_params[0]);
  free(full_params[1]);
  free(full_params[2]);
  free(full_params[3]);
  //Return result of integration
  return 2*PI*result;
}

double forwardIntegrandBinned(double v, void* params)
{
  //Watch out - this doesn't work with any other v-range than 0 to v_max

  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  int N = *((int*)p[1]);
  double v_q = *( (double*)p[2]);
  double* parameters = *((double**)p[3]);

  double beta = v_q/v;

  int bin_number = floor(0.999*N*v/v_max);
  double fv = parameters[bin_number];

  return (1 + asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}


double backwardIntegrandBinned(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  int N = *((int*)p[1]);
  double v_q = *( (double*)p[2]);
  double* parameters = *((double**)p[3]);

  double beta = v_q/v;

  int bin_number = floor(0.999*N*v/v_max);
  double fv = parameters[bin_number];

  return (1 - asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}

double forwardIntegrandPoly(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  int N = *((int*)p[1]);
  double v_q = *( (double*)p[2]);
  double* parameters = *((double**)p[3]);

  double beta = v_q/v;

  double alpha = (v/v_max);

  double logf = 0;

  for(int i = 0; i < N; i++)
  {
    logf -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
  }
  double fv = exp(logf);

  return (1 + asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}

double backwardIntegrandPoly(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  int N = *((int*)p[1]);
  double v_q = *( (double*)p[2]);
  double* parameters = *((double**)p[3]);

  double beta = v_q/v;

  double alpha = (v/v_max);

  double logf = 0;

  for(int i = 0; i < N; i++)
  {
    logf -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
  }
  double fv = exp(logf);

  return (1 - asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}



//Do this more generically???
double multipoleIntegrand(double v, void* params)
{

  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  int N = *((int*)p[1]);
  double v_q = *( (double*)p[2]);
  double* parameters = *((double**)p[3]);


  //std::cout << l << std::endl;
  //std::cout << v_q << std::endl;
  //std::cout << parameters[0] << std::endl;


   double alpha = (v/v_max);

   //double alpha = v/v_max;

  double logf = 0;

  for(int i = 0; i < N; i++)
  {
    //std::cout << parameters[i] << "\t";
    //logf -= pow(alpha,i)*parameters[i];
    logf -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
  }
  //std::cout << std::endl;

  //Check these rates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------------<
  return (v)*exp(logf)*gsl_sf_legendre_Pl (l, v_q/v);

}

/*
double multipoleIntegrand_f0(double v, void* params)
{

  void** p = static_cast<void**>(params);

  //Unpack parameters
  int l = *((int*)p[0]);
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);


  //std::cout << l << std::endl;
  //std::cout << v_q << std::endl;
  //std::cout << parameters[0] << std::endl;


   double alpha = (v/v_max);

   //double alpha = v/v_max;

  double logf1 = 0;
  double logf2 = 0;

  for(int i = 0; i < N_vp; i++)
  {
    //std::cout << parameters[i] << "\t";
    //logf -= pow(alpha,i)*parameters[i];
    logf1 -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
    //logf2 -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i+N_terms];
  }
  //std::cout << std::endl;

  //Check these rates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------------<
  return (v)*(exp(logf1))*gsl_sf_legendre_Pl (l, v_q/v);

}

double multipoleIntegrand_f1(double v, void* params)
{
  void** p = static_cast<void**>(params);

  //Unpack parameters
  int l = *((int*)p[0]);
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);


  //std::cout << l << std::endl;
  //std::cout << v_q << std::endl;
  //std::cout << parameters[0] << std::endl;


   double alpha = (v/v_max);

   //double alpha = v/v_max;

  double logf1 = 0;
  double logf2 = 0;

  int bin = floor(4.999*v/v_max);

  double f2 = parameters[bin+N_terms];

  for(int i = 0; i < N_terms; i++)
  {
    //std::cout << parameters[i] << "\t";
    //logf -= pow(alpha,i)*parameters[i];
    logf1 -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
    //logf2 -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i+N_terms];
  }
  //std::cout << std::endl;

  //Check these rates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------------<
  return (v)*exp(logf1)*(f2)*gsl_sf_legendre_Pl (l, v_q/v);

}
*/

double polyf(double v, void* params)
{
  Astrophysics* astro = ((Astrophysics*)params);
  double alpha = (v/v_max);

  double logf = 0;

  for (int i = 0; i < astro->N_vp; i++)
  {
	  //Change to Chebyshev Polynomials!
   logf -= ChebyshevP(i,2*alpha-1)*astro->vel_params[i];
   //logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*astro->vel_params[i];
  }
  return v*v*exp(logf);
}


//Need to fix this - different normalisation!!!
double polyf_ang(double v, Astrophysics* astro, int k)
{
    double alpha = (v/v_max);

    double logf = 0;

    for (int i = 0; i < astro->N_vp; i++)
    {
  	  //Change to Chebyshev Polynomials!
     logf -= ChebyshevP(i,2*alpha-1)*astro->vel_params_ang[k][i];
     //logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*astro->vel_params[i];
    }
    return exp(logf);	
}

double polyf_total(double v, void* params)
{
    Astrophysics* astro = ((Astrophysics*)params);
    double alpha = (v/v_max);

    double logf = 0;
	double f = 0;
	
	for (int k = 0; k < N_ang; k++)
	{
		logf = 0;
		for (int i = 0; i < astro->N_vp; i++)
		{
		 logf -= ChebyshevP(i,2*alpha-1)*astro->vel_params_ang[k][i];
		}
		f += exp(logf);
	}
    return v*v*f;	
	
	
}

double polyfintegrand_total(double v, void* params)
{
    Astrophysics* astro = ((Astrophysics*)params);
    double alpha = (v/v_max);

    double logf = 0;
	double f = 0;
	
	for (int k = 0; k < N_ang; k++)
	{
		logf = 0;
		for (int i = 0; i < astro->N_vp; i++)
		{
		 logf -= ChebyshevP(i,2*alpha-1)*astro->vel_params_ang[k][i];
		}
		f += exp(logf);
	}
    return v*f;	
	
	
}


/*
double polyf_DIR(double v, void* params)
{

  double alpha = (v/v_max);

  double logf = 0;

  for (int i = 0; i < N_vp; i++)
  {
    //std::cout << ((double*)params)[i] << std::endl;
   //logf -= pow(alpha,i)*((double*)params)[i];
   logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*((double*)params)[i+N_terms];
  }
  return v*v*(exp(logf));
}
*/

double integrator(double v_q, double integrand (double,void*), Astrophysics* astro)
{

  if (v_q > v_max)
  {
    //std::cout << "Hello\t" << v_q << std::endl;
   return 0;
  }

  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

  //Declare gsl function to be integrated
  gsl_function F;
  F.function = integrand;

  F.params = astro;

  double result, error;

  int status = gsl_integration_qag(&F,v_q,v_max, 0, 1e-6, 5000,6,
                             workspace, &result, &error);

  if (status ==  GSL_EROUND)
  {
  //result = 0;
  std::cout << "GSL rounding error!" << std::endl;
  std::cout << result << std::endl;
  }

  //Free workspace
  gsl_integration_workspace_free (workspace);


  //Return result of integration
  return 0.5*4*PI*result;
}


