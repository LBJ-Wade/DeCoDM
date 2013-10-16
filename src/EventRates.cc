
#include <math.h>
#include <iostream>

#include "Detector_Class.h"
#include "ParamSet_Class.h"
#include "EventRates.h"
#include "DMUtils.h"
#include "Parametrisations.h"

#include "gsl/gsl_sf.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"


#ifndef PI
  #define PI 3.14159265358
#endif

//--------File-scope variables------------
double v_max = 1000;


//--------Function Declarations------------
/*
double diffRate (double v, double theta, double phi, ParamSet params)
{
    Detector* expt = params.exptParams;
    double* parameters = params.theoryParams;

    //Experimental parameters
    double m_n = expt->m_n;

    //Theoretical parameters
    double m_x = pow(10,parameters[0]);
    double sigma = pow(10,parameters[1]);
    double v_lagx = parameters[2];
    double v_lagy = parameters[3];
    double v_lagz = parameters[4];
    double sigma_v = parameters[5];

    //Check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   double int_factor = 0;
   if (USE_SD)    int_factor += expt->SD_formfactor(v_min_inverse(abs(v),m_n,m_x))*expt->SD_enhancement();
   if (USE_SI)    int_factor += expt->SI_formfactor(v_min_inverse(abs(v),m_n,m_x))*expt->SI_enhancement();


   return rate_prefactor(m_n, m_x, sigma, 0.3)*int_factor*maxwellRadon(v,theta,phi,v_lagx,v_lagy,v_lagz,sigma_v);
}
*/
double maxwellRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma)
{
  double dotproduct = sin(theta)*(cos(phi)*vlagx + sin(phi)*vlagy) + cos(theta)*vlagz;


 return (pow(2*PI*sigma*sigma,-1.0/2.0))*exp(-0.5*pow((v-dotproduct)/sigma,2));
}

double DMRate(double E, void* params)
{
    Detector* expt = ((ParamSet*)params)->exptParams;
    double* parameters = ((ParamSet*)params)->theoryParams;

    std::vector<double> m_n = expt->m_n;

    double m_x = pow(10,parameters[0]);
    double sigma_SI = pow(10,parameters[1]);
    double sigma_SD = pow(10,parameters[2]);


    //Note the factor of 0.5;

    double rate = 0;

    for (int i = 0; i < expt->N_isotopes; i++)
    {
      double int_factor = 0;
      if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E,i)*expt->SD_enhancement();
      if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E,i)*expt->SI_enhancement(i);

      double v = v_min(E,m_n[i],m_x);

      rate += expt->frac_n[i]*rate_prefactor(m_n[i], m_x, 1, 0.30)*int_factor*currentVelInt(v, params);

    }

    return rate;
}


double spreadRate(double E, void* params)
{
    Detector* expt = ((ParamSet*)params)->exptParams;
    double deltaE = expt->dE;

    //double r = (pow(2*PI*deltaE*deltaE,-0.5))*exp(-0.5*pow((E - currentE),2)/(deltaE*deltaE))*currentRate(E,params);
    double r = (pow(2*PI*deltaE*deltaE,-0.5))*exp(-0.5*pow((E - currentE),2)/(deltaE*deltaE))*currentRate(E,params);
    return r;
}

//--------------------------Preconvolved Rate is currently not being used in the binned likelihood...------------------------

double preConvolvedRate(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
    double deltaE = expt->dE;

    double r = 0.5*(erf((expt->E_max - E)/(sqrt(2)*expt->dE)) - erf((expt->E_min - E)/(sqrt(2)*expt->dE)))*currentRate(E,params);

  //return 0.5*(erf((expt->E_max - E)/(sqrt(2)*expt->dE)) - erf((expt->E_min - E)/(sqrt(2)*expt->dE)))*currentRate(E, params);
  return r;
}

double convolvedRate(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
  if (expt->dE < 1e-3)
  {

    return currentRate(E, params);
  }


  currentE = E;

  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (1000);


  //Declare gsl function to be integrated
  gsl_function F;
  F.function = &spreadRate;

  F.params = params;

  double result, error;

  //switch off default error handler, store old error handler in
  //old_handler
   gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

   double E_lower = currentE - 6*expt->dE;
   double E_higher = currentE + 6*expt->dE;

   if (E_lower < 0) E_lower = 0;
   if (E_higher > 200) E_higher = 200;

    int status = gsl_integration_qag(&F,E_lower,E_higher, 0, 1e-6, 1000,6,
                             workspace, &result, &error);

  if (status ==  GSL_EROUND)
  {
    result = 0;
    std::cout << "GSL rounding error!" << std::endl;
  //std::cout << result << std::endl;
  }

  //Free workspace
  gsl_integration_workspace_free (workspace);

  //Return result of integration
  return result;
}

//-------------------------------------------------------------------------
//-----------Velocity Integrals--------------------------------------------
//-------------------------------------------------------------------------

double VelInt_maxwell(double v, void* params)
{
  double* parameters = ((ParamSet*)params)->theoryParams;
  double v_lag = sqrt(pow(parameters[3],2) + pow(parameters[4],2) + pow(parameters[5],2));
  double v_rms = parameters[6];
  double v_esc = parameters[7];

  double N = 1.0/(pow(2*PI,1.5)*pow(v_rms,3)*gsl_sf_erf(v_esc/(sqrt(2)*v_rms)) - 4*PI*v_rms*v_rms*exp(-0.5*pow(v_esc/v_rms,2)));
  double vel_integral = 0;

  if (v > (v_esc + v_lag))
  {
   return 0;
  }
  else if (v < (v_esc - v_lag))
  {
    vel_integral = sqrt(2)*pow(PI,1.5)*(pow(v_rms,3)/v_lag)*(gsl_sf_erf((v + v_lag)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
    vel_integral -= 4*PI*v_rms*v_rms*exp(-0.5*pow(v_esc/v_rms,2));
  }
  else
  {
    vel_integral = sqrt(2)*pow(PI,1.5)*(pow(v_rms,3)/v_lag)*(gsl_sf_erf((v_esc)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
    vel_integral -= 2*PI*v_rms*v_rms*((v_esc+v_lag - v)/(v_lag))*exp(-0.5*pow(v_esc/v_rms,2));
  }

  //double vel_integral = (1.0/(2*v_lag))*(gsl_sf_erf((v +v_lag)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
  return 2*PI*N*vel_integral;
}

//Returns the unnormalised Lisanti et al directionally averaged speed distribution
double Lisanti_f(double v, void* params)
{
 double v0 = ((double*)params)[0];
 double v_esc = ((double*)params)[1];
 int k = round(((double*)params)[2]);

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
   double v0 = ((double*)params)[0];
  double v_esc = ((double*)params)[1];
  int k = round(((double*)params)[2]);



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
double VelInt_Lisanti(double v, void* params)
{
    double* parameters = ((ParamSet*)params)->theoryParams;

    double velParams[3];
    velParams[0] = parameters[3];
    velParams[1] = parameters[4];
    velParams[2] = parameters[5];

   double vel_integral = integrator(v, &LisantiIntegrand, velParams);


  return vel_integral;

}



//-----------------------------------------------------------------
//-------------Binned parametrisation------------------------------
//-----------------------------------------------------------------


double VelInt_isotropicBinned(double v, void* params )
{
    if (v > v_max) return 0;

    double* parameters = ((ParamSet*)params)->theoryParams;

    double* bin_values;
    bin_values = (double*)malloc(N_terms*sizeof(double));

    double* v_edges;
    v_edges = (double*)malloc((N_terms+1)*sizeof(double));

    //Calculate bin edges
    double v_a = 0;
    double v_b = v_max;
    double delta_v = (v_b - v_a)/N_terms;

    for (int i = 0; i < N_terms+1; i++)
    {
      v_edges[i] = v_a + i*delta_v;
    }

    //Calculate normalised bin values
    for (int i = 0; i < N_terms; i++)
    {
      bin_values[i] = 3*parameters[i+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
      //std::cout << bin_values[i] << "\t";
    }

   double vel_integral = f0_binned(v,bin_values,v_edges, N_terms);

   free(bin_values);
   free(v_edges);

  return vel_integral;
}

double VelInt_forwardBinned(double v, void* params)
{
  double* parameters = ((ParamSet*)params)->theoryParams;

  double* forward_bins;
  forward_bins = (double*)malloc(N_terms*sizeof(double));

  double* backward_bins;
  backward_bins = (double*)malloc(N_terms*sizeof(double));

  double* v_edges;
  v_edges = (double*)malloc((N_terms+1)*sizeof(double));

  //Calculate bin edges
  double v_a = 0;
  double v_b = 1000;

  double delta_v = (v_b - v_a)/N_terms;

    for (int i = 0; i < N_terms+1; i++)
    {
      v_edges[i] = v_a + i*delta_v;
    }

  //Calculate normalised bin values
   for (int i = 0; i < N_terms; i++)
    {
      forward_bins[i] = 3*parameters[i+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
      backward_bins[i] = 3*parameters[i+N_terms+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
    }


    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, forward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, backward_bins);

  free(forward_bins);
  free(backward_bins);
  free(v_edges);

  return vel_integral;
}

double VelInt_backwardBinned(double v, void* params )
{
  double* parameters = ((ParamSet*)params)->theoryParams;

  double* forward_bins;
  forward_bins = (double*)malloc(N_terms*sizeof(double));

  double* backward_bins;
  backward_bins = (double*)malloc(N_terms*sizeof(double));

  double* v_edges;
  v_edges = (double*)malloc((N_terms+1)*sizeof(double));

  //Calculate bin edges
  double v_a = 0;
  double v_b = 1000;

  double delta_v = (v_b - v_a)/N_terms;

    for (int i = 0; i < N_terms+1; i++)
    {
      v_edges[i] = v_a + i*delta_v;
    }

  //Calculate normalised bin values
   for (int i = 0; i < N_terms; i++)
    {
     forward_bins[i] = 3*parameters[i+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
     backward_bins[i] = 3*parameters[i+N_terms+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
    }

    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, backward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, forward_bins);


  free(forward_bins);
  free(backward_bins);
  free(v_edges);

  return vel_integral;
}

//-----------------------------------------------------------------
//-------------Exponential-Polynomial parametrisation--------------
//-----------------------------------------------------------------

double VelInt_isotropicPoly(double v, void* params )
{
    double* parameters = ((ParamSet*)params)->theoryParams;

    //Pack parameters
    double* polyParams;
    polyParams = (double*)malloc(N_terms*sizeof(double));
    for (int i = 0; i < N_terms; i++)
    {
      polyParams[i] = parameters[i+3];
    }

    double vel_integral = multipoleRadon(v, 0, &multipoleIntegrand,polyParams);

    free(polyParams);

    return vel_integral;
}

double VelInt_forwardPoly(double v, void* params )
{
  double* parameters = ((ParamSet*)params)->theoryParams;

  double* forward_bins;
  forward_bins = (double*)malloc(N_terms*sizeof(double));

  double* backward_bins;
  backward_bins = (double*)malloc(N_terms*sizeof(double));

  //Calculate normalised bin values
   for (int i = 0; i < N_terms; i++)
    {
      forward_bins[i] = parameters[i+3];
      backward_bins[i] = parameters[i+N_terms+3];
    }


    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, forward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, backward_bins);

  free(forward_bins);
  free(backward_bins);

  return vel_integral;
}


double VelInt_backwardPoly(double v, void* params )
{
  double* parameters = ((ParamSet*)params)->theoryParams;

  double* forward_bins;
  forward_bins = (double*)malloc(N_terms*sizeof(double));

  double* backward_bins;
  backward_bins = (double*)malloc(N_terms*sizeof(double));

  //Calculate normalised bin values
   for (int i = 0; i < N_terms; i++)
    {
      forward_bins[i] = parameters[i+3];
      backward_bins[i] = parameters[i+N_terms+3];
    }

    //Calculate backward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, backward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, forward_bins);

  free(forward_bins);
  free(backward_bins);

  return vel_integral;
}



double multipoleRadon(double v_q, int l, double integrand (double,void*), double* params)
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

  //Pack paramaters
  void *full_params[3];

  full_params[0] = malloc(sizeof(int));
  full_params[1] = malloc(sizeof(double));

  //Not sure about this bit...

  //-------------------------------------------------------This bit will need editing depending on whether it's 1-D or 3-D------------------------
  full_params[2] = malloc(N_terms*sizeof(double));

  *((int*)full_params[0]) = l;
  *((double*)full_params[1]) = v_q;
  *((double**)full_params[2]) = params;

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

  //Return result of integration
  return 2*PI*result;
}

double forwardIntegrandBinned(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);

  double beta = v_q/v;

  int bin_number = floor(0.999*N_terms*v/v_max);
  double fv = parameters[bin_number];

  return (1 + asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}


double backwardIntegrandBinned(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
    int l = *((int*)p[0]);
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);

  double beta = v_q/v;

  int bin_number = floor(0.999*N_terms*v/v_max);
  double fv = parameters[bin_number];

  return (1 - asin(sqrt(1-beta*beta))*beta/(sqrt(1-beta*beta)))*v*fv;
}

double forwardIntegrandPoly(double v, void* params)
{
  void** p = static_cast<void**>(params);

    //Unpack parameters
  int l = *((int*)p[0]);
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);

  double beta = v_q/v;

  double alpha = (v/v_max);

  double logf = 0;

  for(int i = 0; i < N_terms; i++)
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
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);

  double beta = v_q/v;

  double alpha = (v/v_max);

  double logf = 0;

  for(int i = 0; i < N_terms; i++)
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
  double v_q = *( (double*)p[1]);
  double* parameters = *((double**)p[2]);


  //std::cout << l << std::endl;
  //std::cout << v_q << std::endl;
  //std::cout << parameters[0] << std::endl;


   double alpha = (v/v_max);

   //double alpha = v/v_max;

  double logf = 0;

  for(int i = 0; i < N_terms; i++)
  {
    //std::cout << parameters[i] << "\t";
    //logf -= pow(alpha,i)*parameters[i];
    logf -= gsl_sf_legendre_Pl(i, 2*alpha-1)*parameters[i];
  }
  //std::cout << std::endl;

  //Check these rates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------------<
  return (v)*exp(logf)*gsl_sf_legendre_Pl (l, v_q/v);

}

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

  for(int i = 0; i < N_terms; i++)
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


double polyf(double v, void* params)
{

  double alpha = (v/v_max);

 // double alpha = v/v_max;

  double logf = 0;

  for (int i = 0; i < N_terms; i++)
  {
    //std::cout << ((double*)params)[i] << std::endl;
   //logf -= pow(alpha,i)*((double*)params)[i];
   logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*((double*)params)[i];
  }
  return v*v*exp(logf);
}

double polyf_DIR(double v, void* params)
{

  double alpha = (v/v_max);

 // double alpha = v/v_max;

  double logf = 0;

  for (int i = 0; i < N_terms; i++)
  {
    //std::cout << ((double*)params)[i] << std::endl;
   //logf -= pow(alpha,i)*((double*)params)[i];
   logf -= gsl_sf_legendre_Pl(i,2*alpha-1)*((double*)params)[i+N_terms];
  }
  return v*v*(exp(logf));
}


double integrator(double v_q, double integrand (double,void*), double* params)
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

  F.params = params;

  double result, error;

  int status = gsl_integration_qag(&F,v_q,1000, 0, 1e-6, 5000,6,
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

//------------------------------------------------------------------
//------------Background rates--------------------------------------
//------------------------------------------------------------------

//Background rate in events /keV/kg/day
double BGRate(double E, void* params)
{
 Detector* expt = ((ParamSet*)params)->exptParams;
 double rate = expt->BG_level;

 return rate;

}

