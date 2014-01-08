#include <string.h>
#include <stdio.h>

#include "Detector_Class.h"
#include "ParamSet_Class.h"
#include "EventRates.h"
#include "DMUtils.h"

#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"


#ifndef PI
  #define PI 3.14159265358
#endif

//---------Global variables

std::string expt_folder;
std::string events_folder;
int N_expt;
int mode;
int USE_SD;
int USE_SI;
int dir;
int N_terms;

int USE_FLOAT_BG;
int USE_ASIMOV_DATA;

double (*currentRate) (double, void*);
double (*currentVelInt) (double, void*);
double currentE;

double E_a;
double E_b;

//------------Function Declarations-------------



void setCurrentRate(double rate (double, void*))
{
  currentRate = rate;
}

void setCurrentVelInt(double VelInt (double, void*))
{
 currentVelInt = VelInt;
}

double N_expected(double rate (double,void*), ParamSet parameters)
{
  currentRate = rate;
  //Calculates the total number of expected events
  //using a generic differential event rate ('rate')

  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

	 Detector* expt = parameters.exptParams;
	 double deltaE = expt->dE;

  //switch off default error handler, store old error handler in
  //old_handler
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

  //Declare gsl function to be integrated
  gsl_function F;


  F.params = &parameters;

  double result, error;
  int status;
  if (deltaE < 1e-3)
  {
    F.function = currentRate;

    status = gsl_integration_qag(&F,expt->E_min,expt->E_max, 0.1, 0, 5000,6,
                             workspace, &result, &error);
  }
  else
  {
    F.function = preConvolvedRate;
    E_a = expt->E_min;
    E_b = expt->E_max;
    
    double E_low = E_a - 6*deltaE;
    double E_high = E_b + 6*deltaE;
    if (E_low < 0) E_low = 0;
    if (E_high > 150) E_high = 120;


   status = gsl_integration_qag(&F,E_low,E_high, 0.1, 0, 5000,6,
                             workspace, &result, &error);
  }


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

double N_expected(double rate (double,void*), ParamSet parameters, double E1, double E2)
{
  currentRate = rate;
  //Calculates the total number of expected events
  //using a generic differential event rate ('rate')

  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
         = gsl_integration_workspace_alloc (5000);

	 Detector* expt = parameters.exptParams;
	 double deltaE = expt->dE;

  //switch off default error handler, store old error handler in
  //old_handler
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

  //Declare gsl function to be integrated
  gsl_function F;


  F.params = &parameters;

  double result, error;
  int status;
  if (deltaE < 1e-3)
  {
    F.function = currentRate;

    status = gsl_integration_qag(&F,E1,E2, 0.1, 0, 5000,6,
                             workspace, &result, &error);
  }
  else
  {
    F.function = preConvolvedRate;
    E_a = E1;
    E_b = E2;

    double E_low = E1 - 6*deltaE;
    double E_high = E2 + 6*deltaE;
    if (E_low < 0) E_low = 0;
    if (E_high > 120) E_high = 120;

    status = gsl_integration_qag(&F,E_low,E_high, 0.1, 0, 5000,6,
				 workspace, &result, &error);
  }


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


double rate_prefactor(double m_n, double m_x, double sigma, double rho)
{
  //Calculates the prefactor for the differential event rate for
  //DIRECTIONAL DETECTION of WIMPs, in events day^-1 kg^-1 keV^-1.
  //sigma - interaction cross section in cm^2
  //rho - local DM density in GeV cm^-3
  //m_x - WIMP mass in Gev
  //m_n - nuclear mass in amu

  //Calculate reduced WIMP-proton mass
  double mu = reduced_m(1,m_x);

  //Calculate Prefactor (NB: V0 in km/s)
  double p = 1.274e-12*sigma*rho/(4.0*PI*m_x*mu*mu);

  //Rescale units
  //p *= 1.474e-17*86400;


  return p;
}

double v_min(double E, double m_n, double m_x)
{
 //Calculates the minimum WIMP speed (in km s^-1) required for a recoil of energy E
 //E - recoil energy in keV
 //m_n - nuclear mass in amu
 //m_x - WIMP mass in GeV

 //Calculate reduced mass in GeV
 double mu = reduced_m_GeV(m_n,m_x);
  //double mu = reduced_m(m_n,m_x);

  //double delta = 15;

  //double v = (E*1e-6*(m_n*931.5e-3)/mu + delta*1e-6);
  //v *= 3e5/(sqrt(2*E*1e-6*m_n*931.5e-3));

 //Calculate minimum speed
  return 3e5*sqrt((E/1e6)*(m_n*931.5e-3)/(2*mu*mu));
  //return v;
  //return 5.154e-25*sqrt((E)*(m_n)/(2*mu*mu));

}

double reduced_m(double m_n, double m_x)
{
  //Calculates the reduced mass of the WIMP-nucleus system in kg
  //m_x - WIMP mass in Gev
  //m_n - Nuclear mass in amu

  //Convert nuclear mass into GeV
  m_n *= (931.5e-3);

  //Calculate reduced mass
  double mu = 1.78e-27*(m_n*m_x)/(m_n+m_x);

  //Convert reduced mass to kg
  //mu *= (1e9*1.78e-36);

  return mu;
}

double reduced_m_GeV(double m_n, double m_x)
{
  //Convert nuclear mass into GeV
  m_n *= (931.5e-3);

  //Calculate reduced mass
  double mu = (m_n*m_x)/(m_n+m_x);

  return mu;
}

double p_min(double E, double m_n)
{
  //Calculate p_min(MeV) from the energy E(keV) and nuclear mass m_n(amu)

  return 1e-3*sqrt(m_n*9.315e5*E/2.0);
}

double v_min_inverse(double v, double m_n, double m_x)
{
 //Inverts the calculation of v_min to return the corresponding recoil energy
 //E - recoil energy in keV
 //m_n - nuclear mass in amu
 //m_x - WIMP mass in GeV

 //Calculate reduced mass in
 double mu = reduced_m(m_n,m_x)/(1.78e-27);

 //Calculate corresponding energy
  return 1e6*2*pow(mu*v/3e5,2)/(m_n*9.315e-3);
}

double u_max(double v, double m_n, double m_x)
{
 //Calculates the maximum WIMP speed which can downscatter to below speed v
 //v - speed (in same units as u)
 //m_n - nuclear mass in amu
 //m_x - WIMP mass in GeV

 //Calculate m_n in GeV
  m_n *= 931.5e-3;

 //Calculate maximum speed
  return v*(sqrt(4*m_x*m_n)/(m_x-m_n));
}


double PoissonLike(Detector* expt, ParamSet parameters, double signal_rate (double,void*), double No, double E1, double E2)
{
    //Returns -ln of the poisson likelihood of observing No events when Ne are expected

  //Calculate expected numbers of events
  double Ne = (expt->m_det)*(expt->exposure)*N_expected(signal_rate, parameters, E1, E2);
  double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate,parameters, E1, E2);

  double Ne_tot = Ne+Ne_BG;

  return Ne_tot - No*log(Ne_tot) + logfactNo(No);
}


double signal_fraction(Detector* expt, ParamSet parameters, double signal_rate (double, void*))
{
  //Calculate expected numbers of events
  double Ne = (expt->m_det)*(expt->exposure)*N_expected(signal_rate, parameters);
  double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);

  double Ne_tot = Ne+Ne_BG;

  //Calculate signal fraction
  return Ne/Ne_tot;
}



//double* calcBins(double start, double end, int N_bins)
//{
 //Avoid memory leaks!

//}



double sphericalHarmonic(int l, int m, double theta, double phi)
{

  double azimuthalPart = 0;
  if (m >= 0)
  {
    azimuthalPart = cos(m*phi);
  }
  else
  {
    azimuthalPart = sin(-m*phi);
  }

 return gsl_sf_legendre_sphPlm (l, abs(m), cos(theta))*azimuthalPart;
}

double logfactNo(double No)
{
  /*
 double result = 0;

 for (int i = 0; i < No; i++)
  {
    result += log(i+1);
  }
  return result;
 */
  return 0;
}

int nCr(int n, int r)
{
  //Recursive formula for binomial coefficients
  if (r == 0) return 1;
  if (n == r) return 1;
  return nCr(n-1,r-1) + nCr(n-1,r);
}


//--------------------------------------------------
//---------------File I/O Utilities-----------------
//--------------------------------------------------

int read_param_int(std::ifstream* file, std::string param_name)
{
  std::string name;
  int value;

  if (file->is_open())
  {
    //Reset to start of file
    file->clear();
    file->seekg(0, std::ios::beg);

    while ( (file->good())&&(!file->eof()) )
    {
      *file >> name;


      if (name.compare(param_name)== 0)
      {

	*file >> value;
	return value;
      }
   }
   std::cout << "Parameter not found: " << param_name << std::endl;
   return -1;
  }
  std::cout << "File not open." << std::endl;
  return -1;
}

double read_param_double(std::ifstream* file, std::string param_name)
{
  std::string name;
  double value;

  if (file->is_open())
  {
    //Reset to start of file
    file->clear();
    file->seekg(0, std::ios::beg);

    while ( (file->good())&&(!file->eof()) )
    {
      *file >> name;


      if (name.compare(param_name) == 0)
      {
	*file >> value;
	return value;
      }
   }
   std::cout << "Parameter not found: " << param_name << std::endl;
   return -1;
  }
  std::cout << "File not open." << std::endl;
  return -1;
}

double read_param_vector(std::ifstream* file, std::string param_name, double* output)
{
  std::string name;

  if (file->is_open())
  {
    //Reset to start of file
    file->clear();
    file->seekg(0, std::ios::beg);

    while ( (file->good())&&(!file->eof()) )
    {
      *file >> name;


      if (name.compare(param_name) == 0)
      {
	for (int i = 0; i < 3; i++)
	{
	  *file >> output[i];
	}
	return 1;
      }
   }
   std::cout << "Parameter not found: " << param_name << std::endl;
   return -1;
  }
  std::cout << "File not open." << std::endl;
  return -1;
}

std::string read_param_string(std::ifstream* file, std::string param_name)
{
  std::string name;
  std::string value;


  if (file->is_open())
  {
    //Reset to start of file
    file->clear();
    file->seekg(0, std::ios::beg);

    while ( (file->good())&&(!file->eof()) )
    {
      *file >> name;

      if (name.compare(param_name) == 0)
      {
	*file >> value;
	return value;
      }
   }
   std::cout << "Parameter not found: " << param_name << std::endl;
   return "Error";
  }
  std::cout << "File not open." << std::endl;
  return "Error";
}

int load_params(std::string filename)
{
  //Load global parameters (rho_0 etc.) from filename

  //Open detector parameter file
  std::ifstream file (filename.c_str());
  if (file.is_open())
  {

    //Read in parameter values
    N_expt = read_param_int(&file, "N_expt");
    expt_folder = read_param_string(&file, "expt_folder");
    events_folder = read_param_string(&file, "events_folder");

    mode = read_param_int(&file, "mode");
    dir = read_param_int(&file, "dir");

    USE_SD = read_param_int(&file, "USE_SD");
    USE_SI = read_param_int(&file, "USE_SI");

    USE_FLOAT_BG = read_param_int(&file, "USE_FLOAT_BG");
    USE_ASIMOV_DATA = read_param_int(&file, "USE_ASIMOV_DATA");

    file.close();

    //Display parameter values
    std::cout << "Using parameters:" << std::endl;
    std::cout << "\tN_expt:\t" << N_expt << std::endl;
    std::cout << "\tmode:\t" << mode << std::endl;
    std::cout << "\tdir:\t" << dir << std::endl;
    std::cout << "\tUSE_SI:\t" << USE_SI << std::endl;
    std::cout << "\tUSE_SD:\t" << USE_SD << std::endl;
    std::cout << "\tUSE_FLOAT_BG:\t" << USE_FLOAT_BG << std::endl;
    std::cout << "\tUSE_ASIMOV_DATA:\t" << USE_ASIMOV_DATA << std::endl;
    std::cout <<std::endl;
  }
  else std::cout << "Unable to open experimental parameter file:\t'" << filename << "'" << std::endl;

  return 0;
}


//--------------------------------------
//--------Chebyshev Polynomials---------
//--------------------------------------

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

  return -1;



}
