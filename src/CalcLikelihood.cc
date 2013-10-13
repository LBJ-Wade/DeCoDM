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

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//--------Function Prototypes-----------

double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int mode, int dir);


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
	return 1e30;
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
	  if (experiments[i].USE_BINNED_DATA) experiments[i].bin_data();
	  experiments[i].load_asimov_data(events_folder + "Asimov_Events"+std::string(numstr)+".txt");
    }

  count++;
  }

  theory.m_x = pow(10,params[0])
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
  int offset = USE_SI+USE_SD+N_expt*USE_FLOAT_BG;

  //---------------------------------------------------------------------------------------------------
  //--------Binned Speed parametrisation: mode = 1-------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 1)
  {
      //N_params++;
      double N_vp = N_params-offset;
      astro.initialise_bins(N_vp, DIR);
      astro.vel_params[0] = 0;
      for (int i = 1; i < N_vp; i++)
      {
          astro.vel_params[i] = params[i+offset];
      }
      if (astro.normalise_bins() == -1) return 1e-30;
      astro.rescale_bins(DIR);
      astro.velocityIntegral = &velInt_isotropicBinned;
  }

  //---------------------------------------------------------------------------------------------------
  //--------Binned momentum parametrisation: mode = 2--------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 2)
  {
      std::cout << "Momentum binned is not currently supported!" << std::endl;
      exit (EXIT_FAILURE);
  }
  //---------------------------------------------------------------------------------------------------
  //--------Poly-Exp parametrisation: mode = 3--------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 3)
  {

    //N_params++;
    N_vp = N_params - offset;
    astro.initialise_terms(N_vp, DIR);
    astro.vel_params[0] = 0;
    for (int i = 1; i < N_vp; i++)
    {
          astro.vel_params[i] = params[i+offset];
    }
    astro.normalise_terms(DIR);
    astro.velocityIntegral = &velInt_isotropicPoly;
   }
  

  if (USE_FLOAT_BG)
    {
      for (int i = 0; i < N_expt; i++)
	{
	  experiments[i].BG_level = pow(10,params[1+USE_SD+USE_SI+i]);
	  // std::cout << experiments[i].BG_level << "\t";
	}
      //std::cout << std::endl;
    }

  


  /*
  std::cout << "Input parameters:" << std::endl;
  for (int i = 0; i < *num_hard; i++)
  {
      std::cout << params[i] << std::endl;
  }
  */

  for (int i = 0; i < N_expt; i++)
  {
    loglike += likelihood(&(experiments[i]), &theory, &astro, DIR);
  }

  //    std::cout << -loglike << std::endl << std::endl;

  *result = loglike;
  double LL = 1.0*loglike;
  return LL;
}



double likelihood(Detector* expt , Particlephysics* theory, Astrophysics* astro, int mode, int dir)
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
	    //Calculate expected numbers of events
	    double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
	    double Ne_BG = (expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);

	    double Ne_tot = Ne+Ne_BG;

	    //Calculate signal/BG fractions
	    double f_S = Ne/Ne_tot;
	    double f_BG = Ne_BG/Ne_tot;


	    PL = +Ne_tot - No*log(Ne_tot) + logfactNo(No);


	    //Calculate the event-by-event part
	    double eventLike = 0;
	    for (int i = 0; i < No; i++)
	    {
	      if (dir == 0)  //Do it isotropically!
	      {
		  setCurrentRate(&DMRate);
		  eventLike = f_S*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne;
		  setCurrentRate(&BGRate);
		  eventLike += f_BG*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_BG;
		  PL -= log(eventLike);
	      }
	      else 		//Do it directionally!
	      {
		std::cout << "Warning! Mode 0 does not yet support directional data!" << std::endl;
		//THIS IS DOES NOT WORK
		//PL -= log((expt->exposure)*(expt->m_det)*diffRate(v_min(expt->data[i].energy,expt->m_n,pow(10,((double *) params)[0])),expt->data[i].theta, expt->data[i].phi, parameters)/Ne);
	      }

		  //PL -= log((expt->exposure)*(expt->m_det)*binnedRate(v_min(expt->data[i].energy,expt->m_n,pow(10,((double *) params)[0])),expt->data[i].theta, expt->data[i].phi,5,v_edges,parameters)/Ne);
	    }
	  }

  return PL;

}

