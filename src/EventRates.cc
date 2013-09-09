
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

double maxwellRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma)
{
  double dotproduct = sin(theta)*(cos(phi)*vlagx + sin(phi)*vlagy) + cos(theta)*vlagz;


 return (pow(2*PI*sigma*sigma,-1.0/2.0))*exp(-0.5*pow((v-dotproduct)/sigma,2));
}

double maxwell(double E, void* params)
{
  //double v_max = 1000;

    Detector* expt = ((ParamSet*)params)->exptParams;
    double* parameters = ((ParamSet*)params)->theoryParams;

  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);
  double v_lag = sqrt(pow(parameters[3],2) + pow(parameters[4],2) + pow(parameters[5],2));
  double v_rms = parameters[6];


  double v = v_min(E,m_n,m_x);

  //Note the factor of 0.5;

  double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

  return 0.5*(4*PI)*rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*(1.0/(2*v_lag))*(gsl_sf_erf((v +v_lag)/(sqrt(2)*v_rms)) - gsl_sf_erf((v -v_lag)/(sqrt(2)*v_rms)));
}

double binnedRate(double v, double theta, double phi, int N_bins, double* v_edges, ParamSet params)
{
    Detector* expt = params.exptParams;
    double* parameters = params.theoryParams;

    //Experimental parameters
    double m_n = expt->m_n;

    //Theoretical parameters
    double m_x = pow(10,parameters[0]);
    double sigma_SI = pow(10,parameters[1]);
    double sigma_SD = pow(10,parameters[2]);

    double* bin_values;
    bin_values = (double*)malloc(N_bins*sizeof(double));
    //double* bin_values2;
    //bin_values2 = (double*)malloc(N_bins*sizeof(double));

    for (int i = 0; i < N_bins; i++)
    {
      bin_values[i] = 3*parameters[i+3]/(pow(v_edges[i+1],3) - pow(v_edges[i],3));
    }

    for (int i = 0; i < N_bins; i++)
    {
     // bin_values2[i] = parameters[i+N_bins+2];
    }



    //l = 0
    double vel_integral = 2*PI*(1/((4*PI)))*f0_binned(v,bin_values,v_edges, N_bins);

    //l = 1
    //vel_integral += PI*f1_binned(v,bin_values2,v_edges,N_bins);
    //vel_integral += f1_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(1,-1,theta,phi)*ang_values[1]);
    //vel_integral += f1_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(1,1,theta,phi)*ang_values[2]);

    //l = 2
    //vel_integral += f2_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(2,0,theta,phi)*ang_values[3]);
    //vel_integral += f2_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(2,-2,theta,phi)*ang_values[4]);
    //vel_integral += f2_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(2,-1,theta,phi)*ang_values[5]);
    //vel_integral += f2_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(2,1,theta,phi)*ang_values[6]);
    //vel_integral += f2_binned(v,bin_values,v_edges,N_bins)*(sphericalHarmonic(2,2,theta,phi)*ang_values[7]);

    free(bin_values);
    //free(bin_values2);


   double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(v_min_inverse(abs(v),m_n,m_x))*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(v_min_inverse(abs(v),m_n,m_x))*expt->SI_enhancement();

   return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;


}

//This function has been TAMPERED WITH to work for the momentum method...


double isotropicRateBinned(double E, void* params )
{
    //double v_max = 1000;

    Detector* expt = ((ParamSet*)params)->exptParams;
    double* parameters = ((ParamSet*)params)->theoryParams;



  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

  //std::cout << m_x << "\t" << sigma_SI << "\t" << sigma_SD << std::endl;

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
    //std::cout << std::endl;


  //double v = (reduced_m(40,m_x))*p_min(E,m_n)/(reduced_m(m_n,m_x));
  double v = v_min(E,m_n,m_x);


   double vel_integral = f0_binned(v,bin_values,v_edges, N_terms);


    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

   //std::cout << sigma_SI << std::endl;

   free(bin_values);
   free(v_edges);

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;


}



double spreadRate(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
  double deltaE = expt->dE;


  double r = (pow(2*PI*deltaE*deltaE,-0.5))*exp(-0.5*pow((E - currentE),2)/(deltaE*deltaE))*currentRate(E,params);
  return r;

}

//--------------------------Preconvolved Rate is currently not being used in the binned likelihood...------------------------

double preConvolvedRate(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
/*
  if (expt->dE < 1e-3)
  {
   if ((E > expt->E_min)&&(E < expt->E_max))
   {
    return currentRate(E,params);
   }
   else
   {
    return 0;
   }
  }
  */


  return 0.5*(erf((expt->E_max - E)/(sqrt(2)*expt->dE)) - erf((expt->E_min - E)/(sqrt(2)*expt->dE)))*currentRate(E, params);

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

double forwardRateBinned(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
  double* parameters = ((ParamSet*)params)->theoryParams;

  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

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


  double v = v_min(E,m_n,m_x);

    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, forward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, backward_bins);


    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

  free(forward_bins);
  free(backward_bins);
  free(v_edges);

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;
}





double backwardRateBinned(double E, void* params )
{
  Detector* expt = ((ParamSet*)params)->exptParams;
  double* parameters = ((ParamSet*)params)->theoryParams;

  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

  double* forward_bins;
  forward_bins = (double*)malloc(N_terms*sizeof(double));

  double* backward_bins;
  backward_bins = (double*)malloc(N_terms*sizeof(double));

  double* v_edges;
  v_edges = (double*)malloc((N_terms+1)*sizeof(double));

  //-------------------------------------------------------<<<<<<<<<<<<DONT NEED TO DO THIS FOR ALL BINS - JUST CALCULATE IN THE INTEGRAND
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

  double v = v_min(E,m_n,m_x);

    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandBinned, backward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandBinned, forward_bins);

    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

	 free(forward_bins);
  free(backward_bins);
  free(v_edges);

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;
}



double forwardRatePoly(double E, void* params )
{
  Detector* expt = ((ParamSet*)params)->exptParams;
  double* parameters = ((ParamSet*)params)->theoryParams;

  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

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


  double v = v_min(E,m_n,m_x);

    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, forward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, backward_bins);

	//std::cout << vel_integral << std::endl;

    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

  free(forward_bins);
  free(backward_bins);

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;
}


double backwardRatePoly(double E, void* params )
{
    Detector* expt = ((ParamSet*)params)->exptParams;
  double* parameters = ((ParamSet*)params)->theoryParams;

  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

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


  double v = v_min(E,m_n,m_x);

    //Calculate forward rate
  double vel_integral = multipoleRadon(v, 0, &forwardIntegrandPoly, backward_bins);
  vel_integral += multipoleRadon(v, 0, &backwardIntegrandPoly, forward_bins);

	//std::cout << vel_integral << std::endl;

    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

  free(forward_bins);
  free(backward_bins);

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;
}


double isotropicRatePoly(double E, void* params )
{
    //double v_max = 1000;

    Detector* expt = ((ParamSet*)params)->exptParams;
    double* parameters = ((ParamSet*)params)->theoryParams;



  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);

  double v = v_min(E,m_n,m_x);


  //Note the factor of 0.5;

  //Think about what this factor should actually be...

  // std::cout << N_terms << std::endl;

  //Pack parameters
  double* polyParams;
  polyParams = (double*)malloc(N_terms*sizeof(double));
  for (int i = 0; i < N_terms; i++)
  {
   polyParams[i] = parameters[i+3];
  }



      //l = 0
    double vel_integral = multipoleRadon(v, 0, &multipoleIntegrand,polyParams);

    free(polyParams);

  double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();


  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;



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

//Background rate in events /keV/kg/day
double BGRate(double E, void* params)
{
 Detector* expt = ((ParamSet*)params)->exptParams;
 double rate = expt->BG_level;

 return rate;

}

//Returns the unnormalised Lisanti et al directionally averaged speed distribution
double Lisanti_f(double v, void* params)
{
 double v0 = ((double*)params)[0];
 double v_esc = ((double*)params)[1];
 int k = round(((double*)params)[2]);

 //if (v > 0.5*(sqrt(4*v_esc*v_esc - 3*v0*v0) - v0)) return 0;
 //if (v > (v_esc + v0)) return 0;


 double tot = 0;

 double A = v_esc*v_esc - v*v - v0*v0;
 double B = v0*v;
 double C = k*v0*v0;

 for (int i = 0; i < k; i++)
 {
  tot += (C/B)*pow(-1,i)*(nCr(k, i)/(k-i))*(exp((k-i)*(A+B)/C) - exp((k-i)*(A-B)/C));
  //std::cout << tot << std::endl;
 }
 tot += 2*pow(-1,k);

 if (tot < 0) return 0;

 //std::cout << v0 << "\t" << v_esc << "\t" << k << std::endl;

 //return v*v*pow((exp((v_esc*v_esc - v*v)/(k*v0*v0))-1),k);
 return v*v*tot;
}

double LisantiIntegrand(double v, void* params)
{
  void** p = static_cast<void**>(params);
  double* parameters = *((double**)p[2]);

  double v0 = ((double*)parameters)[0];
  double v_esc = ((double*)parameters)[1];
  int k = round(((double*)parameters)[2]);

  //std::cout << v0 << "\t" << v_esc << "\t" << k << std::endl;




    //if (v > 0.5*(sqrt(4*v_esc*v_esc - 3*v0*v0) - v0)) return 0;
  //if (v > (v_esc + v0)) return 0;


  double tot = 0;

  double A = v_esc*v_esc - v*v - v0*v0;
  double B = v0*v;
  double C = k*v0*v0;

  for (int i = 0; i < k; i++)
  {
    tot += (C/B)*pow(-1,i)*(nCr(k, i)/(k-i))*(exp((k-i)*(A+B)/C) - exp((k-i)*(A-B)/C));
  }
  tot += 2*pow(-1,k);

 if (tot < 0) return 0;
 //std::cout << v0 << "\t" << v_esc << "\t" << k << std::endl;

 //return v*v*pow((exp((v_esc*v_esc - v*v)/(k*v0*v0))-1),k);
 return v*tot/(2*PI);
}

//Returns the norm of the Lisanti distribution (i.e. the integral over all speeds)
double Lisanti_norm(void* params)
{
  double v_esc = ((double*)params)[1];


  //Declare gsl workspace (5000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

  //Declare gsl function to be integrated
  gsl_function F;
  F.function = &Lisanti_f;

  F.params = params;

  double result, error;

    int status = gsl_integration_qag(&F,0,v_max, 0, 1e-6, 5000,6,
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
double LisantiRate(double E, void* params)
{
  Detector* expt = ((ParamSet*)params)->exptParams;
    double* parameters = ((ParamSet*)params)->theoryParams;



  double m_n = expt->m_n;

  double m_x = pow(10,parameters[0]);
  double sigma_SI = pow(10,parameters[1]);
  double sigma_SD = pow(10,parameters[2]);


  double velParams[3];
  velParams[0] = parameters[3];
  velParams[1] = parameters[4];
  velParams[2] = parameters[5];


  //double v = (reduced_m(40,m_x))*p_min(E,m_n)/(reduced_m(m_n,m_x));
  double v = v_min(E,m_n,m_x);


   double vel_integral = multipoleRadon(v, 0, &LisantiIntegrand, velParams);


    double int_factor = 0;
   if (USE_SD)    int_factor += sigma_SD*expt->SD_formfactor(E)*expt->SD_enhancement();
   if (USE_SI)    int_factor += sigma_SI*expt->SI_formfactor(E)*expt->SI_enhancement();

   //std::cout << sigma_SI << std::endl;

  return rate_prefactor(m_n, m_x, 1, 0.3)*int_factor*vel_integral;

}


