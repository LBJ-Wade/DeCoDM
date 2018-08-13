
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
double maxwellRadon(double v, double theta, double phi, Astrophysics* astro)
{
	double total = 0;
    for (int i = 0; i < astro->N_dist; i++)
    {
            double v_rms = astro->v_rms[i];
            double v_esc = astro->v_esc;
            double v_lagx = astro->v_lag_x[i];
			double v_lagy = astro->v_lag_y[i];
			double v_lagz = astro->v_lag_z[i];
	
  		    double dotproduct = sin(theta)*(cos(phi)*v_lagx + sin(phi)*v_lagy) + cos(theta)*v_lagz;

  		  	double N = 1.0/(gsl_sf_erf(v_esc/(sqrt(2)*v_rms)) - sqrt(2.0/PI)*(v_esc/v_rms)*exp(-0.5*pow(v_esc/v_rms,2)));

 		    total += astro->fraction[i]*N*(pow(2*PI*v_rms*v_rms,-1.0/2.0))*(exp(-0.5*pow((v-dotproduct)/v_rms,2)) - exp(-0.5*pow(v_esc/v_rms,2)));
		}
		return total;
}

double maxwellModifiedRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma)
{
  double dotproduct = sin(theta)*(cos(phi)*vlagx + sin(phi)*vlagy) + cos(theta)*vlagz;
  double vlagsq = vlagx*vlagx + vlagy*vlagy + vlagz*vlagz;

 return (2.0*sigma*sigma + vlagsq - dotproduct*dotproduct)*(pow(2*PI*sigma*sigma,-1.0/2.0))*exp(-0.5*pow((v-dotproduct)/sigma,2))/(3e5*3e5);
}

double DMRateNR(double E, void* params)
{
	//Currently just calculates the separate contibutions F_ij^NN', with no overall scaling...
    Detector* expt = ((ParamSet*)params)->exptParams;
    Particlephysics* theory = ((ParamSet*)params)->theoryParams;
    Astrophysics* astro = ((ParamSet*)params)->astroParams;

    std::vector<double> m_n = expt->m_n;

    //Note the factor of 0.5;

	int op1 = theory->op1;
	int op2 = theory->op2;
	
	int N1 = theory->N1;
	int N2 = theory->N2;

    double rate1 = 0;
	double total_rate = 0;

	double A = 0;
	double B = 0;

	 double conv = pow((1.98e-14*1.0/(theory->m_x+0.9315)),2)/(16.0*PI);
	  
	  int i_c = 0;
	  if ((N1 == 0) and (N2 == 0)) i_c = 0;
	  if ((N1 == 0) and (N2 == 1)) i_c = 1;
	  if ((N1 == 1) and (N2 == 0)) i_c = 2;
	  if ((N1 == 1) and (N2 == 1)) i_c = 3;

    for (int i = 0; i < expt->N_isotopes; i++)
    {
      double int_factor = 0;
      //Define conversion factor from amu-->keV
      double amu = 931.5*1000;
      //Convert recoil energy to momentum transfer q in keV
      double  q = sqrt(2*m_n[i]*amu*E);
      //Convert q into fm^-1
      //q *= (1e-12/1.97e-7);

      double v = v_min(E,m_n[i],theory->m_x);

	  rate1 = 0;
	  A = 0;
	  B = 0;
	  
	  double Eta = astro->velocityIntegral(v, astro);
	  double MEta = astro->modifiedVelocityIntegral(v, astro);
	  
	  //Currently using the Fitzpatrick et al normalisation...
	  
	  if (op1 == op2)
	  {
		  switch (op1) 
		  {
			  case 1:
			      rate1 = expt->F_M(E, i, i_c)*Eta;
				  break;
			  case 3:
			  	  A = MEta*expt->F_Sigma1(E, i, i_c);
			      B = 0.25*pow((q/amu),2.0)*Eta*expt->F_Phi2(E, i, i_c);
			      rate1 = pow((q*1e-6),2.0)*(A+B);
				  break;
			  case 4:
			  	  //Use the correct values here...Check the normalisation/scaling...
			      //rate1 = (1.0/16.0)*Eta*(expt->F_Sigma1(E, i, i_c) + expt->F_Sigma2(E, i, i_c));
			      rate1 = (1.0/4.0)*(1.0/16.0)*Eta*expt->F_SD(E, i, i_c);
				  break;
			  case 5:
			  	  A = MEta*expt->F_M(E, i, i_c);
			      B = Eta*pow((q/amu),2.0)*expt->F_Delta(E, i, i_c);
			      rate1 = pow((q*1e-6),2.0)*(A + B)/4.0;
				  break;
			  case 6:
				  rate1 = pow((q*1e-6),4.0)*Eta*expt->F_Sigma2(E, i, i_c)/16.0;
				  break;
			  case 7:
				  rate1 = MEta*expt->F_Sigma1(E, i, i_c)/8.0;
				  break;
			  case 8:
				  A = MEta*expt->F_M(E, i, i_c);
				  B = Eta*pow((q/amu),2.0)*expt->F_Delta(E, i, i_c);
				  rate1 = (A + B)/4.0;
				  break;
			  case 9:
				  rate1 = Eta*pow((q*1e-6),2.0)*expt->F_Sigma1(E, i, i_c)/16.0;
				  break;
			  case 10:
				  rate1 = Eta*pow((q*1e-6),2.0)*expt->F_Sigma2(E, i, i_c)/4.0;
				  break;
			  case 11:
				  rate1 = Eta*pow((q*1e-6),2.0)*expt->F_M(E, i, i_c)/4.0;
				  break;
				  
			  //Long range interactions
			  case 101:
			      rate1 = pow(q*1e-6,-4.0)*expt->F_M(E, i, i_c)*Eta;
				  break;
			  case 104:
				  rate1 = 0;
				  break;
			  case 105:
				  A = MEta*expt->F_M(E, i, i_c);
			      B = Eta*pow((q/amu),2.0)*expt->F_Delta(E, i, i_c);
			      rate1 = pow((q*1e-6),-2.0)*(A + B)/4.0;
				  break;
			   case 106:
				  rate1 = Eta*expt->F_Sigma2(E, i, i_c)/16.0;
				  break;
   			   case 111:
			      rate1 = Eta*pow((q*1e-6),-2.0)*expt->F_M(E, i, i_c)/4.0;
   				  break;  
				  
		  }
	  }
	  else
	  {
		  if ((op1 == 8 && op2 == 9)||(op1 == 9 && op2 == 8))
		  {
		  	rate1 = pow((q*1e-6),2.0)*(1.0/(amu*1e-6))*Eta*expt->F_Sigma1Delta(E, i, i_c)/8.0;	
		  }
		  //if (op1 < 100) return 0;
		  //if (op2 < 100) return 0;
		  //These are long range only
		  if ((op1 == 104 && op2 == 105)||(op1 == 105 && op2 == 104))
		  {
			  rate1 = -(1.0/(amu*1e-6))*Eta*expt->F_Sigma1Delta(E, i, i_c)/8.0;
		  }
		  if ((op1 == 104 && op2 == 106)||(op1 == 106 && op2 == 104))
		  {
			  rate1 = Eta*expt->F_Sigma2(E, i, i_c)/16.0;
		  }
		  
		  //std::cout << "No interference terms have been included yet!" << std::endl;
	  }
	  total_rate += expt->frac_n[i]*rate1*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*conv;
	  
    }
	return expt->efficiency(E)*total_rate;
}

double DMRate(double E, void* params)
{
    Detector* expt = ((ParamSet*)params)->exptParams;
    Particlephysics* theory = ((ParamSet*)params)->theoryParams;
    Astrophysics* astro = ((ParamSet*)params)->astroParams;

    std::vector<double> m_n = expt->m_n;

	//theory->PrintAll();

    //Note the factor of 0.5;
    double rate = 0;

    for (int i = 0; i < expt->N_isotopes; i++)
    {
      double int_factor = 0;

      //Currently only isoscalar SD scattering (i.e. 0 rather than 1 in the 'component' field)
      if ((USE_SD)&&(pow(expt->J[i],2) > 1e-6))	int_factor += theory->sigma_SD*(16.0/3.0)*(1.0/16.0)*(expt->F_Sigma1(E, i, 0)+expt->F_Sigma2(E, i, 0));
	  //if (USE_SI)    int_factor += theory->sigma_SI*expt->SI_formfactor(E,i)*expt->SI_enhancement(i);
	  if (USE_SI)    
	  {
		  double sig = 0;
		  sig = (1.973e-14*1.973e-14)*4.0*pow(reduced_m_GeV(1.0, theory->m_x),2.0)/PI;

		  int_factor += 0.5*sig*expt->SI_formfactor(E,i)*
		  (pow(theory->lambda_p_D*expt->N_p[i] + theory->lambda_n_D*expt->N_n[i],2.0)
			  + pow(theory->lambda_p_Dbar*expt->N_p[i] + theory->lambda_n_Dbar*expt->N_n[i],2.0));
	  }
	 
      //Define conversion factor from amu-->keV
      double amu = 931.5*1000;
      //Convert recoil energy to momentum transfer q in keV
      double  q = sqrt(2*m_n[i]*amu*E);
      //Convert q into fm^-1
      q *= (1e-12/1.97e-7);

      double v = v_min(E,m_n[i],theory->m_x);
	  double R1 = expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*int_factor*astro->velocityIntegral(v, astro);
	  //Add in the O7 operator - need to check/include interaction factor...
	  /*
	  double p5 = 0;
	  double R5 = 0;
	  if (theory->sigma_O5 > 1e-60)
	  {
	      p5 = theory->sigma_O5*(expt->F_M(E, i, 0)*astro->modifiedVelocityIntegral(v, astro) + pow(q/amu, 2.0)*expt->F_Delta(E,i,0)*astro->velocityIntegral(v, astro));
	      R5 = (1.0/4.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*(pow(q/amu, 2.0))*p5;
	  }
	  
	  double R2 = 0;
	  if (theory->sigma_O7 > 1e-60)
	  {
	      R2 =  (1.0/8.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*theory->sigma_O7*expt->F_Sigma1(E, i, 0)*astro->modifiedVelocityIntegral(v, astro);
	  }
	  
	  double p3 = 0;
	  double R3 = 0;
	  if (theory->sigma_O15 > 1e-60)
	  {
	      p3 = theory->sigma_O15*(expt->F_Sigma1(E, i, 0)*astro->modifiedVelocityIntegral(v, astro) + 2*pow(q/amu, 2.0)*expt->F_Phi2(E,i,0)*astro->velocityIntegral(v, astro));
	      R3 = (1.0/32.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*(pow(q/amu, 4.0))*p3;
      }
	  
	  rate += R1 + R2 + R3 + R5;
	  */
	  rate += R1;
    }
    return rate;
}

double DMRateDirectional(double E, double theta, double phi, void* params)
{
    Detector* expt = ((ParamSet*)params)->exptParams;
    Particlephysics* theory = ((ParamSet*)params)->theoryParams;
    Astrophysics* astro = ((ParamSet*)params)->astroParams;

    std::vector<double> m_n = expt->m_n;

    //Note the factor of 0.5;

    //theory->PrintAll();

    double rate = 0;

    for (int i = 0; i < expt->N_isotopes; i++)
    {
      double int_factor = 0;

	  //Do a list of nuclear response functions for all the operators...
	  //Also do a list of proton correction factors for all operators...

      //Currently only isoscalar SD scattering (i.e. 0 rather than 1 in the 'component' field)

      if ((USE_SD)&&(pow(expt->J[i],2) > 1e-6))	int_factor += theory->sigma_SD*(16.0/3.0)*(1.0/16.0)*(expt->F_Sigma1(E, i, 0)+expt->F_Sigma2(E, i, 0));
	  //if (USE_SI)    int_factor += theory->sigma_SI*expt->SI_formfactor(E,i)*expt->SI_enhancement(i);
	  if (USE_SI)    
	  {
		  double sig = 0;
		  sig = (1.973e-14*1.973e-14)*4.0*pow(reduced_m_GeV(1.0, theory->m_x),2.0)/PI;
		  int_factor += 0.5*sig*expt->SI_formfactor(E,i)*
		  (pow(theory->lambda_p_D*expt->N_p[i] + theory->lambda_n_D*expt->N_n[i],2.0)
			  + pow(theory->lambda_p_Dbar*expt->N_p[i] + theory->lambda_n_Dbar*expt->N_n[i],2.0));
	  }

      //Define conversion factor from amu-->keV
      double amu = 931.5*1000;
      //Convert recoil energy to momentum transfer q in keV
      double  q = sqrt(2*m_n[i]*amu*E);
      //Convert q into fm^-1
      q *= (1e-12/1.97e-7);

      double v = v_min(E,m_n[i],theory->m_x);
	  //std::cout << std::endl;
	  double R1 = 0;
	  for (int i = 0; i < astro->N_dist; i++)
	  {
		  R1 += astro->fraction[i]*maxwellRadon(v, theta, phi, astro);
	  }
      R1 *= expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*int_factor;
	  //std::cout << R1 << std::endl;
	  //Need to generalise to any velocity distribution...
	  
	  double p5 = 0;
	  double R5 = 0;
	  if (theory->sigma_O5 > 1e-60)
	  {
	      p5 = theory->sigma_O5*(expt->F_M(E, i, 0)*maxwellModifiedRadon(v, theta, phi, 0, 0, 220.0, 156.0) + pow(q/amu, 2.0)*expt->F_Delta(E,i,0)*maxwellRadon(v, theta, phi, astro));
	      R5 = (1.0/4.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*(pow(q/amu, 2.0))*p5;
	  }
	  
	  double R2 = 0;
	  if (theory->sigma_O7 > 1e-60)
	  {
	      R2 = (1.0/8.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*theory->sigma_O7*expt->F_Sigma1(E, i, 0)*maxwellModifiedRadon(v, theta, phi, 0, 0, 220.0, 156.0);
	  }
	  //std::cout << R1 << "\t" << R2 << std::endl;
	  
	  double p3 = 0;
	  double R3 = 0;
	  if (theory->sigma_O15 > 1e-60)
	  {
	      p3 = theory->sigma_O15*(expt->F_Sigma1(E, i, 0)*maxwellModifiedRadon(v, theta, phi, 0, 0, 220.0, 156.0) + 2*pow(q/amu, 2.0)*expt->F_Phi2(E,i,0)*maxwellRadon(v, theta, phi, astro));
	      R3 = (1.0/32.0)*expt->frac_n[i]*rate_prefactor(m_n[i], theory->m_x, 1, Astrophysics::rho_x)*(pow(q/amu, 4.0))*p3;
	  }
	  rate += R1 + R2 + R3 + R5;
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

