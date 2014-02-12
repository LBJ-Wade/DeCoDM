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

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//--------Function Prototypes-----------

double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int vmode, int dir);


//------------Interface with fortran--------

extern "C" {
double loglikelihood_(double * params, int* num_hard, double *result);
}


//------------Function definitions----------


double loglikelihood_(double * params, int* num_hard, double *result)
{
  double loglike = 0;

  static int count = 0;

  static std::vector<Detector> experiments;
  Astrophysics astro; //NEED TO WORRY ABOUT STATIC-NESS
  Particlephysics theory;

  //Initialise and load in data on first run
  if (count == 0)
  {
    //Load in general parameters
    load_params("params.ini");
    astro.load_params();

    //Check to make sure SI and/or SD inteactions are being used
    if ((USE_SI + USE_SD) < 1)
    {
	std::cout << "Must specify SI and/or SD interactions in params.ini" << std::endl;
	exit (EXIT_FAILURE);
    }

    //Check to make sure correct 'direction' option is being used
    if ((DIR != 0)&&(DIR != 1))
    {
	std::cout << "Must specify dir = 0 or dir = 1 in params.ini" << std::endl;
        exit (EXIT_FAILURE);
    }


    if ((DIR == 1))
    {
	std::cout << "Directional detection is not currently supported..." << std::endl;
        exit (EXIT_FAILURE);
    }

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
    }

  count++;
  }

  //Load in 'theory' parameters from input sample
  theory.m_x = pow(10,params[0]);
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


  int N_params = *num_hard;
  int offset = USE_SI+USE_SD+N_expt*USE_FLOAT_BG+ 6*USE_VARY_FF;

  //---------------------------------------------------------------------------------------------------
  //--------Benchmark speed distributions: vmode = 0----------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (vmode == 0)
  {
    //Use parameters from 'dist.txt' file
    astro.load_params();
  }

  //---------------------------------------------------------------------------------------------------
  //--------Binned Speed parametrisation: vmode = 1----------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  else if (vmode == 1)
  {
      double N_vp = N_params-offset;
      astro.initialise_bins(N_vp, DIR);
      astro.vel_params[0] = 0;
      for (int i = 1; i < N_vp; i++)
      {
          astro.vel_params[i] = params[i+offset];
      }

      //Exit if bins cannot be normalised
      if (astro.normalise_bins() == -1) return 1e-30;

      astro.rescale_bins(DIR);
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
  //---------------------------------------------------------------------------------------------------
  //--------Poly-Exp parametrisation: vmode = 3--------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  else if (vmode == 3)
  {

    //N_params++;
    int N_vp = N_params - offset - 1;
    astro.initialise_terms(N_vp, DIR);
    //astro.vel_params[0] = 0;
    for (int i = 0; i < N_vp; i++)
    {
          astro.vel_params[i] = params[i+offset+1];
    }
    //astro.normalise_terms(DIR);
    astro.velocityIntegral = &velInt_isotropicPoly;
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
  }


  //for (int i = 0; i < experiments[0].alpha.size(); i++)
  //{
    //std::cout << experiments[0].N[0] << "\t" << experiments[0].N[2] << std::endl;
    //std::cout << experiments[0].alpha[0] << "\t" << experiments[0].alpha[2] << std::endl;
    //std::cout << experiments[0].beta[0] << "\t" << experiments[0].beta[2] << std::endl;
  //}
  //std::cout << std::endl;
  /*
  std::cout << "Input parameters:" << std::endl;
  for (int i = 0; i < *num_hard; i++)
  {
      std::cout << params[i] << std::endl;
  }
  */

  for (int i = 0; i < N_expt; i++)
  {
    loglike += likelihood(&(experiments[i]), &theory, &astro, vmode, DIR);
  }

  *result = loglike;
  double LL = 1.0*loglike;
  return LL;
}



double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int vmode, int dir)
{

    ParamSet parameters(expt,theory, astro);

    int No = expt->data.size();
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
	    //std::cout << "rho_x:\t" << astro->rho_x << std::endl;

	    //Calculate expected numbers of events
	    double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
	    double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);

	    double Ne_tot = Ne+Ne_BG;

	    //Calculate signal/BG fractions
	    double f_S = Ne/Ne_tot;
	    double f_BG = Ne_BG/Ne_tot;



	    PL = +Ne_tot - No*log(Ne_tot) + logfactNo(No);

	    //std::cout << "Ne_S:\t" << Ne << ";\tNe_BG:\t" << Ne_BG << std::endl;


	    //Calculate the event-by-event part
	    double eventLike = 0;
	    for (int i = 0; i < No; i++)
	    {
	      if (dir == 0)  //Do it isotropically!
	      {
		  setCurrentRate(&DMRate);
		  eventLike = (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
		  setCurrentRate(&BGRate);
		  eventLike += (expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_tot;
		  PL -= log(eventLike);
	      }
	      else 		//Do it directionally!
	      {

	      }
	    }
	  }

  return PL;

}

