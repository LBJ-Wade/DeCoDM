
#include <math.h>
#include <iostream>

#include "Detector_Class.h"
#include "ParamSet_Class.h"
#include "EventRates.h"
#include "DMUtils.h"


#include "Astrophysics_Class.h"
#include "Particlephysics_Class.h"

#include "gsl/gsl_sf.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"


#ifndef PI
  #define PI 3.14159265358
#endif

//--------File-scope variables------------


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
    Particlephysics* theory = ((ParamSet*)params)->theoryParams;
    Astrophysics* astro = ((ParamSet*)params)->astroParams;

    std::vector<double> m_n = expt->m_n;

    //Note the factor of 0.5;

    double rate = 0;

    for (int i = 0; i < expt->N_isotopes; i++)
    {
      double int_factor = 0;

      //Currently only isoscalar SD scattering (i.e. 0 rather than 1 in the 'component' field)
      if ((USE_SD)&&(pow(expt->J[i],2) > 1e-6))	int_factor += theory->sigma_SD*expt->SD_formfactor(E,i, 1)*expt->SD_enhancement(i);
      if (USE_SI)    int_factor += theory->sigma_SI*expt->SI_formfactor(E,i)*expt->SI_enhancement(i);

      double v = v_min(E,m_n[i],theory->m_x);

      rate += expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*int_factor*astro->velocityIntegral(v, astro);

    }
    return rate;
}
//Need to work out how best to do this in GSL!! Nested integrals suck!
/*
double DMRate_timedep(double E, void* params)
{
    Detector* expt = ((ParamSet*)params)->exptParams;
    Particlephysics* theory = ((ParamSet*)params)->theoryParams;
    Astrophysics* astro = ((ParamSet*)params)->astroParams;
    
    std::vector<double> m_n = expt->m_n;
    
    //Note the factor of 0.5;
    
    double rate = 0;
    
    for (int i = 0; i < expt->N_isotopes; i++)
    {
        double int_factor = 0;
        
        //Currently only isoscalar SD scattering (i.e. 0 rather than 1 in the 'component' field)
        if ((USE_SD)&&(pow(expt->J[i],2) > 1e-6))	int_factor += theory->sigma_SD*expt->SD_formfactor(E,i, 1)*expt->SD_enhancement(i);
        if (USE_SI)    int_factor += theory->sigma_SI*expt->SI_formfactor(E,i)*expt->SI_enhancement(i);
        
        double v = v_min(E,m_n[i],theory->m_x);
        
        double modulation = (1 + astro->mod_amplitude*cos(2.0*PI*(time - astro->mod_phase)/astro->mod_period))
        
        rate += expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*int_factor*astro->velocityIntegral(v, astro);
        
    }
    
    
    return rate;
}
*/
 
 

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

    double r = 0.5*(erf((E_b - E)/(sqrt(2)*expt->dE)) - erf((E_a - E)/(sqrt(2)*expt->dE)))*currentRate(E,params);

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

