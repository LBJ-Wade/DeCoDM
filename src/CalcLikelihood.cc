#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>

#include "Event_Class.h"
#include "Detector_Class.h"
#include "ParamSet_Class.h"
#include "Astrophysics_Class.h"
#include "Particlephysics_Class.h"
#include "EventRates.h"
#include "DMUtils.h"
#include "Distributions.h"
#include "Neutrinos.h"
#include "Parametrisations.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


#ifndef PI
  #define PI 3.14159265358
#endif

//--------Function Prototypes-----------

double loglikelihood(double * params, int* num_hard, double *result);
double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int vmode, int dir);
double geteventnumbers(double * params, int* num_hard, double *Ne);
double RecalcAsimovData(double m_x, double sigma_SI, double sigma_SD);
double RecalcData(double m_x, double sigma_SI, double sigma_SD);
double CalcAsimovLike(double m_x, double sigma_SI, double sigma_SD,double A_day, double A_night);
double CalcDirLike(double m_x, double sigma_SI, double sigma_SD);
double CalcMixedLike(double m_x, double Ne, double A);
void generateEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD);

//static std::vector<Detector> experiments;

//------------Interface with fortran--------

extern "C" double loglikelihood_(double * params, int* num_hard, double *result)
{
	return loglikelihood(params, num_hard, result);
}

extern "C" double geteventnumbers_(double * params, int* num_hard, double *Ne)
{
    return geteventnumbers(params, num_hard, Ne);
}

extern "C" { void dsinterface_nevents2_( double*, double*, double*, double*, double*, double*, double* ); }

//------------Interface with python---------


//------------Function definitions----------


double geteventnumbers(double * params, int* num_hard, double *Ne)
{
	
	
	
	
}

//-----------------------------------------------------------------------------



double loglikelihood(double * params, int* num_hard, double *result)
{
  double loglike = 0;

  static int count = 0;

  
  static std::vector<Detector> experiments;



  Astrophysics astro; //NEED TO WORRY ABOUT STATIC-NESS
  Particlephysics theory;

  //Initialise and load in data on first run
  if (count == 0)
  {
	  //Initialise neutrino tables
	  //LoadFluxTable();
	  
    //Load in general parameters
    load_params("params.ini");
    astro.load_params();

    //Check to make sure SI and/or SD inteactions are being used
    if ((USE_SI + USE_SD) < 1)
    {
	std::cout << "Must specify SI and/or SD interactions in params.ini" << std::endl;
	exit (EXIT_FAILURE);
    }

	/*
    if ((DIR == 1))
    {
		std::cout << "Directional detection is not currently supported..." << std::endl;
		exit (EXIT_FAILURE);
    }*/

    char numstr[21]; // enough to hold all numbers up to 64-bits

      //Load experimental parameters
    for (int i = 0; i < N_expt; i++)
    {
	  sprintf(numstr, "%d", i+1);
	  experiments.push_back(Detector(expt_folder + "Experiment"+std::string(numstr)+".txt"));
	  experiments[i].load_data(events_folder + "Events"+std::string(numstr)+".txt");
	  if (experiments[i].USE_BINNED_DATA)
	  {
	    experiments[i].bin_data();
	    experiments[i].load_asimov_data(events_folder + "Asimov_Events"+std::string(numstr)+".txt");
	  }
	  if (experiments[i].USE_DIR == 1) experiments[i].angular_bin_data(N_ang);
    }

    count++;
  }

  /*
  static int runcount = 0;
  std::cout << runcount << std::endl;
  //if (runcount == 1) exit(0);
  runcount++;	
  */

  //Load in 'theory' parameters from input sample
  theory.m_x = pow(10,params[0]);
  //std::cout << theory.m_x << std::endl;
  if (USE_SI)
  {
    theory.sigma_SI = pow(10,params[1]);

    if (USE_SD)
    {
      theory.sigma_SD = pow(10,params[2]);
    }
    else
    {
      theory.sigma_SD = 0;
    }
  }
  else
  {
    theory.sigma_SI = 0;
    theory.sigma_SD = pow(10,params[1]);
  }
  
  //theory.PrintAll();


  int N_params = *num_hard;
  int offset = USE_SI+USE_SD+N_expt*USE_FLOAT_BG;

  if (N_expt == 1) offset += 6*USE_VARY_FF;
  if (N_expt == 2) offset += 6*USE_VARY_FF;
  if (N_expt == 3) offset += 9*USE_VARY_FF;

  //---------------------------------------------------------------------------------------------------
  //--------Benchmark speed distributions: vmode = 0----------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (vmode == 0)
  {
    //Use parameters from 'dist.txt' file
    astro.load_params();
	if (N_params > 2)
	{
		astro.v_lag[0] = params[2];
		astro.v_lag_z[0] = params[2];
		astro.v_rms[0] = params[3];
	}

  }

  //---------------------------------------------------------------------------------------------------
  //--------Binned Speed parametrisation: vmode = 1----------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  else if ((vmode == 1))
  {
    //std::cout << N_params << std::endl;
      double N_vp = N_params-offset;
      astro.initialise_bins(N_vp, 0);
      astro.vel_params[0] = 0;
      for (int i = 1; i < N_vp; i++)
      {
          astro.vel_params[i] = params[i+offset];
	  //std::cout << astro.vel_params[i] << std::endl;

      }

      //Exit if bins cannot be normalised
      if (astro.normalise_bins() == -1) return 1e30;

      astro.rescale_bins(0);
      astro.velocityIntegral = &velInt_isotropicBinned;
  }

  //---------------------------------------------------------------------------------------------------
  //--------Binned momentum parametrisation: vmode = 2--------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  else if (vmode == 2)
  {

      std::cout << "Momentum binned is not currently supported!" << std::endl;
      exit (EXIT_FAILURE);
  }
  
  //---------------------------------------------------------------------
  //--------Poly-Exp parametrisation: vmode = 3--------------------------      
  //---------------------------------------------------------------------
  else if ((vmode == 3))
  {

    //N_params++;
    //int N_vp = N_params - offset - 1;
	int N_vp = N_terms;
	astro.initialise_terms(N_vp, 1);
    //astro.vel_params[0] = 0;
	//astro.vel_params_ang[0][0] = 0;
	/*
	for (int j = 1; j < N_vp; j++)
	{
			astro.vel_params_ang[0][j] = params[j + offset];
		
		}
	*/
	for (int k = 0; k < N_ang; k++)
	{
		astro.vel_params_ang[k][0] = 0;
	    double a0 = 0;
		for (int j = 1; j < N_vp; j++)
	    {
			  //astro.vel_params_ang[k][j] = params[j + k*(N_vp) + offset];
			astro.vel_params_ang[k][j] = params[j + k*(N_vp-1) + offset];
			a0 -= (astro.vel_params_ang[k][j])*pow(-1.0, j);
	    }
		astro.vel_params_ang[k][0] = a0;
	}
	
	/*
	for (int k = 0; k < N_ang; k++)
	{
		std::cout << " k = " << k << std::endl;
		for (int i = 0; i < N_vp; i++)
		{
			std::cout << "     i = " << i << ": " << astro.vel_params_ang[k][i] << std::endl;		
		}	
	}
	std::cout << std::endl;
	*/
	
	//HERE - sort out normalisation!
	
	double *angnorms = new double[N_ang];
	
    astro.normalise_terms(1, angnorms);
	
	for (int k = 0; k < N_ang; k++)
	{
		params[offset + (N_vp-1)*N_ang + k + 1] = angnorms[k];
	}
	
	delete[] angnorms;
	
	//std::cout << astro.vel_params_ang[0][0] << std::endl;
	//std::cout << astro.vel_params_ang[1][0] << std::endl;
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

	//F.params = &params;
    F.function = &polyf_angint;

	double v1, v2;
	
	
	int N_bins = 100;
	double dv = 1000.0/N_bins;
	astro.N_vp = N_bins;
	
	double *vpars = new double[N_terms];
	
	std::vector<double> binvals;
	astro.bin_params.clear();
	for (int k = 0; k < N_ang; k++)
	{
		binvals.clear();
		for (int j = 0; j < N_terms; j++)
		{
			vpars[j] = astro.vel_params_ang[k][j];
		}
		for (int vi = 0; vi < N_bins; vi++)
		{
			result = 0;
			v1 = vi*dv;
			v2 = (vi+1)*dv;
			F.params = vpars;
	    	status = gsl_integration_qag(&F,v1,v2, 0, 1e-6, 5000,4,
					workspace, &result, &error);
			if (status ==  GSL_EROUND)
			{
				std::cout << " Rounding error in bin calculation..." << std::endl;
				
			}
			binvals.push_back(result/dv);
		}
		astro.bin_params.push_back(binvals);
	}
	
	
	
	delete[] vpars;
	astro.initialise_bins(N_bins, 0);
	
    //astro.velocityIntegral = &velInt_polytotal;
	astro.velocityIntegral = &velInt_dirBinned;
	
	double sumtotal = 0;
	
		/*
		std::cout << " Values of bin parameters..." << std::endl;
		for (int k = 0; k < N_ang; k++)
		{
			std::cout << " k = " << k << std::endl;
			for (int i = 0; i < N_bins; i++)
			{
				sumtotal += astro.bin_params[k][i];
				std::cout << "\t" << astro.bin_params[k][i] << std::endl;
			}
		
		}
		*/

   }
  
  else if ((vmode == 4))
    {
      int N_vp = 1000;
      astro.initialise_bins(N_vp, 0);
      //astro.vel_params[0] = 0;
      for (int i = 0; i < N_vp; i++)
	{
          astro.vel_params[i] = params[i+offset+1];
	  //std::cout << astro.vel_params[i] << std::endl;
	}
      //std::cout << std::endl;

      //Exit if bins cannot be normalised
      //if (astro.normalise_bins() == -1) return 1e30;

      //astro.rescale_bins(DIR);
      astro.velocityIntegral = &velInt_isotropicBinned;
    }
  

  //Load in BG values for each experiment
  if (USE_FLOAT_BG)
    {
      for (int i = 0; i < N_expt; i++)
	{
	  experiments[i].BG_level = pow(10,params[1+USE_SD+USE_SI+i]);
	  //std::cout << experiments[i].BG_level << "\t";
	}
      //std::cout << std::endl;
    }


  //Load in Form Factor parameters for each experiment - only using S_00 at the moment...
  if (USE_VARY_FF)
  {
    for (int k = 0; k < 2; k++)
    {

     experiments[0].N[2*k] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*k + 1];
     experiments[0].alpha[2*k] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*k + 2];
     experiments[0].beta[2*k] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*k + 3];
    }

    experiments[2].N[0] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*2 + 1];
    experiments[2].alpha[0] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*2 + 2];
    experiments[2].beta[0] = params[USE_SD+USE_SI + N_expt*USE_FLOAT_BG + 3*2 + 3];


  }

  /*
  //for (int i = 0; i < experiments[0].alpha.size(); i++)
    //{
    std::cout << experiments[0].N[0] << "\t" << experiments[0].N[2] << std::endl;
    std::cout << experiments[0].alpha[0] << "\t" << experiments[0].alpha[2] << std::endl;
    std::cout << experiments[0].beta[0] << "\t" << experiments[0].beta[2] << std::endl;
    //}
  std::cout << std::endl;

  std::cout << experiments[2].N[0] << std::endl;
  std::cout << experiments[2].alpha[0] << std::endl;
  std::cout << experiments[2].beta[0] << std::endl;

  std::cout << std::endl;
  */



  //Write a file explaining what inputs were used

  /*
  std::ofstream file("f.txt");
  if (file.is_open())
    {  
       std::cout << *num_hard << std::endl;
       std::cout << "Input parameters:" << std::endl;
       for (int i = 0; i < astro.N_vp; i++)
       {
           std::cout << params[i] << std::endl;
           file << astro.vel_params[i] << std::endl;
       }
       std::cout << std::endl;
       file.close();
       exit (EXIT_FAILURE);       
    }
  else std::cout << "Unable to open file f.txt" << std::endl;
  */


  for (int i = 0; i < N_expt; i++)
  {
    loglike -= likelihood(&(experiments[i]), &theory, &astro, vmode, experiments[i].USE_DIR);
	//std::cout << theory.m_x << "\t" << theory.sigma_SI << "\t" << loglike << std::endl;
  }

  
  //This is where I'm quiting MARK123
  std::cout << " Likelihood = " << loglike << std::endl;
  exit(0);
  *result = loglike;
  double LL = 1.0*loglike;
  return LL;
}

//Maybe this wants to be in the CalcLikelihood file - it definitely does

double RecalcAsimovData(double m_x, double sigma_SI, double sigma_SD)
{
	static int count = 0;
	//Check to make sure the detector data has been loaded from file...
	if (count == 0)
	{
        //Load global parameters
        load_params("params.ini");
	
	    char numstr[21]; // enough to hold all numbers up to 64-bits

	      //Load experimental parameters
	    for (int i = 0; i < N_expt; i++)
	    {
		  sprintf(numstr, "%d", i+1);
		  experiments.push_back(Detector(expt_folder + "Experiment"+std::string(numstr)+".txt"));
	    }
		
	    count++;
	}
	
	double A_daynight = 0.00;
	//double A_night = 0.96;
	
	
    Astrophysics astro;

    astro.load_params();
	
    Particlephysics theory;
	Particlephysics underlying;
    theory.m_x = m_x;
    theory.sigma_SI = sigma_SI;
    theory.sigma_SD = sigma_SD;
	
	underlying.m_x = 100;
	underlying.sigma_SI = 1e-46;
	underlying.sigma_SD = 0;

    double scaling = 1;
    if (astro.dist_type == "lisanti")
      {
         scaling = 1.0/(Lisanti_norm(&astro));
      }
	
	for (int j = 0; j < N_expt; j++)
	{
		Detector* expt = &(experiments[j]); 
		ParamSet parameters(expt,&theory, &astro);
		ParamSet param_underlying(expt, &underlying, &astro);
	    //Calculate ASIMOV data if the bin width is defined
	    if ((expt->bin_width > 1e-3))
	    {
	        double Ne_bin = 0;

	        for (int i = 0; i < expt->N_Ebins; i++)
		    {
			  //First need to reset!!!
			  expt->asimov_data[i] = 0;
			  Ne_bin = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters,expt->bin_edges[i], expt->bin_edges[i+1]));
			  expt->asimov_data[i] += Ne_bin;
			  //Ne_bin = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, param_underlying,expt->bin_edges[i], expt->bin_edges[i+1]));
			  //expt->asimov_data[i] += Ne_bin;
			  //std::cout << Ne_bin << std::endl;
		      Ne_bin = expt->m_det*expt->exposure*(N_expected(&BGRate, parameters, expt->bin_edges[i], expt->bin_edges[i+1]));
		      expt->asimov_data[i] += Ne_bin;
			  //expt->asimov_data[i] += expt->neutrino_data[i]*(1 + pow(-1,j)*A_daynight);
			  
			  //std::cout << expt->neutrino_data[i]*(1 + pow(-1,j)*A_daynight) << std::endl;
			}
			
			//std::cout << expt->asimov_data[0] << std::endl;
			//generateBGEvents_Asimov(expt);
			//generateNeutrinoEvents_Asimov(expt);
	    }
		
	}
	
}

double RecalcData(double m_x, double sigma_SI, double sigma_SD)
{
	static int count = 0;
	//Check to make sure the detector data has been loaded from file...
	if (count == 0)
	{
        //Load global parameters
        load_params("params.ini");
	
	    char numstr[21]; // enough to hold all numbers up to 64-bits

	      //Load experimental parameters
	    for (int i = 0; i < N_expt; i++)
	    {
		  sprintf(numstr, "%d", i+1);
		  experiments.push_back(Detector(expt_folder + "Experiment"+std::string(numstr)+".txt"));
	    }
		
	    count++;
	}
	
	for (int i = 0; i < N_expt; i++)
	{
		generateEvents(&(experiments[i]), m_x, sigma_SI,sigma_SD);
	}
	
}


double CalcAsimovLike(double m_x, double sigma_SI, double sigma_SD, double A_day, double A_night)
{
	//NOTE - I need to maximise the two likelihoods separately!!!
	
	//Note that the neutrino fluxes are calculated only when the file is opened,
	//so this is fine - it only affects the likelihood
	//experiments[0].A_nu = A_day*(1+A_night);
	experiments[0].A_nu = A_day;
	//experiments[1].A_nu = A_day*(1-A_night);
	//if (experiments.size() > 1)  experiments[1].A_nu = A_flux*(1.0-A_daynight);

	double loglike = 0;
	//Be careful, I'm doing something a little funny here!
	
	double A_day_true = 1.00;
	double A_night_true = 0.00;

    Astrophysics astro;

    astro.load_params();
	
    Particlephysics theory;
    theory.m_x = m_x;
    theory.sigma_SI = sigma_SI;
    theory.sigma_SD = sigma_SD;

    //double scaling = 1;
    //if (astro.dist_type == "lisanti")
    //  {
    //     scaling = 1.0/(Lisanti_norm(&astro));
    //  }
    for (int i = 0; i < N_expt; i++)
    {
      loglike += likelihood(&(experiments[i]), &theory, &astro, 0, 0);
    }
	//Now include likelihood for the neutrino fluxes
	double sigma_flux = 0.61/5.69;
	double sigma_day = sqrt(sigma_flux*sigma_flux + 0.33*0.33);
	double sigma_night = sigma_day;
	//-----------------------------NOT SURE WHICH VALUES TO USE---
	//-----Roughly!!! - could be as small as 2%, not 5%
	sigma_day = 1.0*sigma_flux;
	sigma_day = 1.0;
	sigma_night = 0.01;
	//sigma_day = 0.001;
	//sigma_night = sigma_day;
	//sigma_night = sqrt(2)*0.61/5.69;
	
	//----Use 16% uncertainty
	//sigma_day = 0.16;
	
	//-------------------------------FIX THIS NUMBER!!!
	//double sigma_daynight = 0.01; //1% uncertainty
	//std::cout << sigma_day << std::endl;
	//std::cout << A_day << "\t" << A_day_true << "\t" << sigma_day << std::endl;
	double L_day = (pow(A_day - A_day_true,2.0))/(2.0*sigma_day*sigma_day);
	double L_night = (pow(A_night - A_night_true,2.0))/(2.0*sigma_night*sigma_night);
	loglike += L_day;
	if (A_day < 0) loglike += 1e30;
	//loglike += L_night;
	//theory.sigma_SI = 1e-60;
	//theory.sigma_SD = 1e-60;
	
	//double L0 = 0;
	//Be careful, I'm doing something a little funny here!
    //for (int i = 0; i < N_expt; i++)
    //{
    //  L0 += likelihood(&(experiments[i]), &theory, &astro, 0, 0);
    //}

    //return 2*(L0-loglike);
	return loglike;
	
}

double CalcDirLike(double m_x, double sigma_SI, double sigma_SD)
{	
    
	double loglike = 0;
	//Be careful, I'm doing something a little funny here!

    Astrophysics astro;

    astro.load_params();
	
	//astro.velocityIntegral = &velInt_maxwell;
	
    Particlephysics theory;
    theory.m_x = m_x;
    theory.sigma_SI = sigma_SI;
    theory.sigma_SD = sigma_SD;

    //double scaling = 1;
    //if (astro.dist_type == "lisanti")
    //  {
    //     scaling = 1.0/(Lisanti_norm(&astro));
    //  }
	
	if (INCLUDE_NU) LoadFluxTable();
	
    for (int i = 0; i < N_expt; i++)
    {
      loglike += likelihood(&(experiments[i]), &theory, &astro, 0, 1);
    }
	if (INCLUDE_NU) ClearFluxTable();
	
	return loglike;
}

double CalcMixedLike(double m_x, double Ne, double A)
{	
    
	if (Ne < 0) return 1e30;
	//if (Ne < 1e-45) return 1e30;
	//if (A < 1e-36) return 1e30;
	//if (A < 0) return 1e30;
	//if (Ne < 1e-6) Ne = 0;
	//if (A < 1e-3) A = 0;
	//if (Ne < 1e-5) Ne = 0;
	if ((A < 0)||(A > 1)) return 1e30;
	//if ((A < 0)) return 1e30;
	//if (A < 1e-5) A = 0;
	//if ((m_x < 49)||(m_x > 51)) return 1e30;
	if (m_x > 1e5) return 1e30;
	//if (A < 1e-6) A = 0;
	//std::cout << A << std::endl;
	double loglike = 0;
	//Be careful, I'm doing something a little funny here!

    Ne = experiments[0].No();

    Astrophysics astro;

    astro.load_params();
	
	//astro.velocityIntegral = &velInt_maxwell;
	
    Particlephysics theory_SI;
    theory_SI.m_x = m_x;
    theory_SI.sigma_SI = 1e-45;
    theory_SI.sigma_SD = 0;
	theory_SI.sigma_O5 = 0;
	theory_SI.sigma_O7 = 0;
	theory_SI.sigma_O15 = 0;
	
    Particlephysics theory_O7;
    theory_O7.m_x = m_x;
    theory_O7.sigma_SI = 0;
    theory_O7.sigma_SD = 0;
	theory_O7.sigma_O5 = 1e-45;
	theory_O7.sigma_O7 = 0;	
	theory_O7.sigma_O15 = 0;
	
	Detector* expt = &(experiments[0]);
	
	ParamSet parameters_SI(expt, &theory_SI, &astro);
    ParamSet parameters_O7(expt, &theory_O7, &astro);

	double Ne_SI = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters_SI);
    double Ne_O7 = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters_O7);
	//std::cout << Ne_SI << "\t" << Ne_O7 << std::endl;
	
	//if (m_x > 1e5)
	//{
	    //std::cout << m_x << "\t" << Ne << "\t" << A <<  std::endl;
	//}
	//std::cout << Ne_SI << "\t" << Ne_O7 << "\t" << Ne << "\t" << A <<  std::endl;
	

	Particlephysics theory_test;
	theory_test.m_x = m_x;
	theory_test.sigma_SI = (Ne_O7/Ne_SI)*A;
	theory_test.sigma_SD = 0;
	theory_test.sigma_O5 = 0;
	theory_test.sigma_O7 = 0;
	theory_test.sigma_O15 = 0;
	
	Particlephysics theory;
	theory.m_x = m_x;
	//theory.sigma_SI = 1e-45*(Ne/Ne_SI)*(1-A);
	theory.sigma_SI = 1e-45*(Ne/Ne_SI)*(1-A);
	if (A > 1) 
	{
		//std::cout << "SI\t" << theory.sigma_SI << std::endl;
		//std::cout << "O7\t" << theory.sigma_O7 << std::endl;
	}
	//theory.sigma_SI = 0;
    theory.sigma_SD = 0;
	//theory.sigma_O15 = 0;
	theory.sigma_O7 = 0;
	//theory.sigma_O15 = 1e-45*(Ne/Ne_O7)*A;
	//theory.sigma_O7 = 1e-40*(Ne/Ne_O7)*(A);
	theory.sigma_O5 = 1e-45*(Ne/Ne_O7)*A;
	//std::cout << theory.sigma_O7 << std::endl;
	//theory.sigma_O7 = 0;
	
	if (Ne <  1e-44)
	{
		std::cout << m_x << "\t" << "Ne_SI: " << Ne << " ; Ne_O7: " << A << std::endl;
		std::cout << likelihood(&(experiments[0]), &theory, &astro, 0, 1) << "\t" << likelihood(&(experiments[0]), &theory_test, &astro, 0, 1) << std::endl;
	}
	
    //double scaling = 1;
    //if (astro.dist_type == "lisanti")
    //  {
    //     scaling = 1.0/(Lisanti_norm(&astro));
    //  }
	
	if (INCLUDE_NU) LoadFluxTable();
	
    for (int i = 0; i < N_expt; i++)
    {
      loglike += likelihood(&(experiments[i]), &theory, &astro, 0, experiments[i].USE_DIR);
    }
	if (INCLUDE_NU) ClearFluxTable();
	
	//if (m_x > 1e3) loglike+=(m_x - 1e3);
	
	return loglike;
}







double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int vmode, int dir)
{

    ParamSet parameters(expt,theory, astro);

    
      //Calculate number of expected and observed events

    double PL = 0;

    if (expt->USE_BINNED_DATA)
    {
	   for (int i = 0; i < expt->N_Ebins; i++)
	   {
	     double No_bin = 0;

	     if (USE_ASIMOV_DATA)
	     {
		No_bin = expt->asimov_data[i];
	     }
	     else
	     {
		No_bin = expt->binned_data[i];
	     }
		 PL += PoissonLike(expt, parameters, &DMRate, No_bin, expt->bin_edges[i], expt->bin_edges[i+1]);

	   }
	}
	else
	{
		if ((vmode != 3)||(dir == 0))
		{
			int No = expt->data.size();
			//Also - maybe I'll use the approximate velocity integral here too!
			//Reset the proper velocity integral...
			if (vmode == 3) 
			{
				//astro->velocityIntegral = &velInt_polytotal;
				astro->velocityIntegral = &velInt_dirBinned;
				parameters.astroParams = astro;
			}
			
		    //Calculate expected numbers of events
		    double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
		    double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);
	        double Ne_nu = 1e-10;
			if (INCLUDE_NU) Ne_nu = (expt->m_det)*(expt->exposure)*N_expected(&NeutrinoRate, parameters);
		
		    double Ne_tot = Ne+Ne_BG+Ne_nu;
	        //std::cout << Ne_tot << std::endl;
		    //Calculate signal/BG fractions
		    double f_S = Ne/Ne_tot;
		    double f_BG = Ne_BG/Ne_tot;

			//std::cout <<  theory->m_x << "\t" << theory->sigma_SI << "\t" << Ne << "\t" << No << "\t" << std::endl;

		    PL = +Ne_tot - No*log(Ne_tot) + logfactNo(No);
			//std::cout << Ne_tot << std::endl;
		    //Calculate the event-by-event part
		    double eventLike = 0;
		    for (int i = 0; i < No; i++)
		    {
				//std::cout << expt->data[i].energy << "\t" << expt->data[i].theta << "\t" <<  expt->data[i].phi << std::endl;
				eventLike = 0;
		      if (dir == 0)  //Do it isotropically!
		      {
				  //---------FIX THE CONVOLVED RATE---------------
			  
				  setCurrentRate(&DMRate);
				  eventLike = (expt->exposure)*(expt->m_det)*DMRate(expt->data[i].energy,&parameters)/Ne_tot;
			  
				  setCurrentRate(&BGRate);
				  eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
				  if (INCLUDE_NU)
				  {
				      setCurrentRate(&NeutrinoRate);
				      eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
			      }
				  PL -= log(eventLike);
				  //std::cout << PL << std::endl;
		      }
		      else 		//Do it directionally!
		      {
				  //Note - currently can't convolve...
				  //setCurrentRate(&DMRate);
				  eventLike = (expt->exposure)*(expt->m_det)*DMRateDirectional(expt->data[i].energy,expt->data[i].theta, expt->data[i].phi,&parameters)/Ne_tot;
			      //std::cout << eventLike << std::endl;
				  //Assuming all backgrounds are isotropic
				  setCurrentRate(&BGRate);
				  eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
				  if (INCLUDE_NU)
				  {
				      setCurrentRate(&NeutrinoRate);
				      eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
			      }
				  PL -= log(eventLike);
		      }
		    }
		}
		else if ((vmode == 3)&&(dir == 1))
		{
			PL = 0;
			for (int j = 0; j < N_ang; j++)
			{
				//std::cout << j << std::endl;
				int No = expt->data_ang[j].size();
			
			    astro->velocityIntegral = &velInt_DRT_disc;
				//astro->velocityIntegral = &velInt_DRT;
				j_bin = j;
				parameters.astroParams = astro;
			    //Calculate expected numbers of events
			    double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
		
		
				std::cout << " j = " << j << ": No = " << No << " ; Ne = " << Ne << std::endl;
		
				
			    char numstr[3]; // enough to hold all numbers up to 64-bits
				sprintf(numstr, "%d", j);
		
				//-----Print out velocity integral
				std::ofstream outfile;
				std::string filename("IRT"+std::string(numstr)+".txt");
				outfile.open(filename.c_str());
			
				for (int iv = 0; iv < 1000; iv++)
				{
					double v = iv*1.0;
					outfile << v << "\t" << astro->velocityIntegral(v, astro) << "\n";
				}	
				outfile.close();
				
		
			    double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);
		        double Ne_nu = 1e-10;
				if (INCLUDE_NU) Ne_nu = (expt->m_det)*(expt->exposure)*N_expected(&NeutrinoRate, parameters);
		
			    double Ne_tot = Ne+Ne_BG+Ne_nu;
		        //std::cout << Ne_tot << std::endl;
			    //Calculate signal/BG fractions
			    double f_S = Ne/Ne_tot;
			    double f_BG = Ne_BG/Ne_tot;

				//std::cout <<  theory->m_x << "\t" << theory->sigma_SI << "\t" << Ne << "\t" << No << "\t" << std::endl;

			    PL += +Ne_tot - No*log(Ne_tot) + logfactNo(No);
				//std::cout << Ne_tot << std::endl;
			    //Calculate the event-by-event part
			    double eventLike = 0;
			    for (int i = 0; i < No; i++)
			    {
					  eventLike = 0;

					  setCurrentRate(&DMRate);
					  eventLike = (expt->exposure)*(expt->m_det)*DMRate(expt->data_ang[j][i].energy,&parameters)/Ne_tot;
		  
					  setCurrentRate(&BGRate);
					  eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data_ang[j][i].energy,&parameters)/Ne_tot;
					  if (INCLUDE_NU)
					  {
					      setCurrentRate(&NeutrinoRate);
					      eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data_ang[j][i].energy,&parameters)/Ne_tot;
				      }
					  PL -= log(eventLike);
			      }
		  	}
			//std::cout << std::endl;	
		}
	  }
  return PL;

}

