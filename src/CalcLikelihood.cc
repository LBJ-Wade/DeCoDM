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
#include "EventRates.h"
#include "DMUtils.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//--------Function Prototypes-----------

double likelihood(Detector* expt , double* params, int mode, int dir);


//------------Interface with fortran--------

extern "C" {
double loglikelihood_(double * params, int* num_hard, double *result);
}


//------------Function definitions----------


double loglikelihood_(double * params, int* num_hard, double *result)
{
  double loglike = 0;

  static int count = 0;

  double norm = 0;
  double norm2 = 0;

  //int N_bins = 0;
  N_terms = 0;

  static std::vector<Detector> experiments;



  //Initialise and load in data on first run
  if (count == 0)
  {
    //Load in general parameters
    load_params("params.ini");

    //Check to make sure SI and/or SD inteactions are being used
    if ((USE_SI + USE_SD) < 1)
    {
	std::cout << "Must specify SI and/or SD interactions in params.ini" << std::endl;
	return 1e30;
    }

    //Check to make sure correct 'direction' option is being used
    if ((dir != 0)&&(dir != 1))
    {
	std::cout << "Must specify dir = 0 or dir = 1 in params.ini" << std::endl;
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

  int N_params = *num_hard;
  int offset = USE_SI+USE_SD+N_expt*USE_FLOAT_BG;

  //---------------------------------------------------------------------------------------------------
  //--------Binned Speed parametrisation: mode = 1-------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 1)
  {
      N_params++;

      N_terms = N_params-1-offset;

      for (int i = 1+offset; i < (N_terms)+offset; i++)
	{
	  //	  std::cout << params[i] << std::endl;
	  //      std::cout << params[i + (N_bins/2) -1] << std::endl;
	  norm += params[i];
	  //norm2 += params[i+(N_bins/2)-1];
	}
  }
  //---------------------------------------------------------------------------------------------------
  //--------Binned momentum parametrisation: mode = 2--------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 2)
  {


      N_terms = N_params-2;

      N_params++;


      for (int i = 2; i < N_terms+1; i++)
      {
	norm += params[i];
      }
      std::cout << "Momentum binned is not currently supported!" << std::endl;
  }
  //---------------------------------------------------------------------------------------------------
  //--------Poly-Exp parametrisation: mode = 3--------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  if (mode == 3)
  {

    N_params++;

    if (dir == 0)
    {
	//Declare gsl workspace (1000 subintervals)
	gsl_integration_workspace * workspace
	      = gsl_integration_workspace_alloc (3000);


	  //Declare gsl function to be integrated
	  gsl_function F;

	  //----------------------------------<<<<<<<<<<<<<<<<<<<FIX OFFSET!!!!!!!!!!!!!!_______________-------------!

	  N_terms = N_params-1-offset;

	  //std::cout << N_terms << std::endl;
	  double* parameters;
	  parameters = (double*)malloc((N_terms)*sizeof(double));

	  for (int i = 1; i < (N_terms); i++)
	  {
	    parameters[i] = params[i+offset];
	  }

	  parameters[0] = 0;

	  F.function = &polyf;
	  F.params = parameters;

	  double error;
	  int status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,workspace, &norm, &error);

	  if (status ==  GSL_EROUND)
	  {
	    std::cout << "GSL rounding error!" << std::endl;
	    std::cout << norm << std::endl;
	  }

      //Free workspace
      gsl_integration_workspace_free (workspace);
      free(parameters);

    }
    else if (dir == 1)
    {
	  N_terms= N_params-1-offset;

	  //Declare gsl workspace (1000 subintervals)
	  gsl_integration_workspace * workspace
	      = gsl_integration_workspace_alloc (3000);


	  //Declare gsl function to be integrated
	  gsl_function F;


	    //std::cout << N_terms << std::endl;
	  double* parameters1;
	  parameters1 = (double*)malloc((N_terms/2)*sizeof(double));

	  double* parameters2;
	  parameters2 = (double*)malloc((N_terms/2)*sizeof(double));

	  parameters1[0] = 0;
	  parameters2[0] = 0;

	  for (int i = 1; i < (N_terms/2); i++)
	  {
	    parameters1[i] = params[i+offset];
	    parameters2[i] = params[i+offset +(N_terms/2)-1];
	    //std::cout << parameters1[i] << "\t" << parameters2[i] << std::endl;
	  }

	  double error;
	  int status;

	  F.function = &polyf;
	  F.params = parameters1;

	  status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,workspace, &norm, &error);

	  if (status ==  GSL_EROUND)
	  {
	    std::cout << "GSL rounding error!" << std::endl;
	    std::cout << norm << std::endl;
	  }

	  F.function = &polyf;
	  F.params = parameters2;

	  status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,workspace, &norm2, &error);

	  if (status ==  GSL_EROUND)
	  {
	    std::cout << "GSL rounding error!" << std::endl;
	    std::cout << norm2 << std::endl;
	  }

	  //Free workspace
	  gsl_integration_workspace_free (workspace);
	  free(parameters1);
	  free(parameters2);
    }
  }

  //Pack up all the input parameters
  double* full_params;
  int N_fullparams = N_params+(1-USE_SI)+(1-USE_SD) - N_expt*USE_FLOAT_BG; //add in extra parameters for un-used interactions
  full_params = (double*)malloc(N_fullparams*sizeof(double));

  //Pack up m_x, sigma_si, sigma_sd
  full_params[0] = params[0];

  if (USE_SI)
  {
    full_params[1] = params[1];

    if (USE_SD)
    {
      full_params[2] = params[2];
    }
    else
    {
      full_params[2] = -100;
    }
  }
  else
  {
    full_params[1] = -100;
    full_params[2] = params[1];
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

  if (mode == 1)
  {
    //std::cout << norm << "\t" << norm2 << std::endl;
    if (norm < 1)
    {
      full_params[3] = 1 - norm;
    }
    else
    {
      *result = 1e30;
      return 1e30;
    }


      for (int i = 4; i < (N_fullparams); i++)
	  {
	    full_params[i] = params[i-3+offset];
	  }

  }
  else if (mode == 3)
  {
    if (dir == 0)
    {
      full_params[3] = log(norm);
      for (int i = 4; i < (N_fullparams); i++)
      {
	full_params[i] = params[i-3+offset];
      }
    }
    else
    {
      full_params[3] = log(norm) - log(params[*num_hard-1]);
      full_params[3+(N_terms/2)] = log(norm2) - log(0.5-params[*num_hard-1]);

      for (int i = 4; i < (N_fullparams-N_terms/2); i++)//---------------------------<<<<<<<<<<<<<NEED TO CHECK ALL THE 'OFFSET' TERMS
      {
	full_params[i] = params[i-3+offset];
	full_params[i+N_terms/2] = params[i+ N_terms/2-4+offset];
      }
    }
  }


  /*
  std::cout << "Input parameters:" << std::endl;
  for (int i = 0; i < *num_hard; i++)
  {
      std::cout << params[i] << std::endl;
  }
  */

  /*
  std::cout << "Full paramset:" << std::endl;
  for (int i = 0; i < 3; i++)
  {
  std::cout << full_params[i] << "\t";
  }
  std::cout << std::endl;

  for (int i = 3; i < 8; i++)
  {
  std::cout << full_params[i] << "\t";
  }
  std::cout << std::endl;
  */

  /*
  for (int i = 8; i < 13; i++)
  {
  std::cout << full_params[i] << "\t";
  }
  std::cout << std::endl <<std::endl;
  */


  for (int i = 0; i < N_expt; i++)
  {
    loglike += likelihood(&(experiments[i]), full_params, mode, dir);
  }

  free(full_params);

  //    std::cout << -loglike << std::endl << std::endl;

  *result = loglike;
  double LL = 1.0*loglike;
  return LL;
}



double likelihood(Detector* expt, double* params, int mode, int dir)
{

    ParamSet parameters(expt,params);

    int No = expt->data.size();
      //Calculate number of expected and observed events

    double PL = 0;

    //---------------------------------------------------------------------------------------------------
    //--------Exact velocity distribution known: mode = 0------------------------------------------------
    //---------------------------------------------------------------------------------------------------
    if (mode == 0)
    {
	int N_dist;
	double fraction;
	double sigma_v;
	double v_lag[3];
	char numstr[21]; // enough to hold all numbers up to 64-bits
	//--------------------------------------------------<<<<<<<<<<<<<<<<<<CURRENTLY ONLY WORKS FOR A SINGLE DISTRIBUTION N_dist=1

	//Read in velocity distribution from dist.txt
	//Open file for reading in distribution parameters


	std::ifstream file ("dist.txt");
	if (file.is_open())
	{
	  //Read in parameter values

	  N_dist = read_param_int(&file, "N_dist");

	  for (int i = 0; i < N_dist; i++)
	  {
	    sprintf(numstr, "%d", i+1);
	    fraction = read_param_double(&file, "fraction"+std::string(numstr));
	    read_param_vector(&file, "v_lag"+std::string(numstr),v_lag);
	    sigma_v = read_param_double(&file, "sigma_v" + std::string(numstr));
	  }


	    // sigma_v = 156;
	  // double v_lag = 230;

	  double full_params[7];
	  full_params[0] = params[0];
	  full_params[1] = params[1];
	  full_params[2] = params[2];
	  full_params[3] = v_lag[0];
	  full_params[4] = v_lag[1];
	  full_params[5] = v_lag[2];
	  full_params[6] = sigma_v;

	  parameters.theoryParams = full_params;

	  setCurrentVelInt(&VelInt_maxwell);

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

	  }
    	else std::cout << "Unable to open distribution parameter file:\t'" << "dist.txt" << "'" << std::endl;

    }

    //---------------------------------------------------------------------------------------------------
    //--------Binned Speed parametrisation: mode = 1-----------------------------------------------------
    //---------------------------------------------------------------------------------------------------
    if (mode == 1)
    {
      if (dir == 0)
      {
	setCurrentVelInt(&VelInt_isotropicBinned);
	//Calculate expected numbers of events...
	double Ne = (expt->m_det)*(expt->exposure)*(N_expected(&DMRate, parameters));
	double Ne_BG = (expt->m_det)*(expt->exposure)*(N_expected(&BGRate, parameters));
	double Ne_tot = Ne+Ne_BG;

	//Calculate signal/BG fractions
	double f_S = Ne/Ne_tot;
	double f_BG = Ne_BG/Ne_tot;

	//Calculate poisson part of the log-likelihood
	PL = +Ne_tot - No*log(Ne_tot) + logfactNo(No);

	 if (isnan(PL))
	  {
	    return 1e30;
	  }

	//Calculate the event-by-event part
	double eventLike = 0;
	for (int i = 0; i < No; i++)
	  {
	    setCurrentRate(&DMRate);
	    eventLike = f_S*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne;
	    setCurrentRate(&BGRate);
	    eventLike += f_BG*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_BG;
	    PL -= log(eventLike);
	  }
	  if (isnan(PL))
	  {
	    return 1e30;
	  }
      }
      else if (dir == 1)
      {

	//Calculate expected number of events
	setCurrentVelInt(&VelInt_forwardBinned);
	double Ne_forward = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
	setCurrentVelInt(&VelInt_backwardBinned);
	double Ne_backward = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);

	//Expected background assuming isotropy
	double Ne_BG_forward = 0.5*(expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);
	double Ne_BG_backward = Ne_BG_forward;

	double Ne_tot_forward = Ne_forward + Ne_BG_forward;
	double Ne_tot_backward = Ne_backward + Ne_BG_backward;

	//Calculate signal/BG fractions
	double f_S_forward = Ne_forward/(Ne_tot_forward);
	double f_BG_forward = 1 - f_S_forward;
	double f_S_backward = Ne_backward/(Ne_tot_backward);
	double f_BG_backward = 1 - f_S_backward;


      int No_forward = 0;
      int No_backward = 0;

      for (int i = 0; i < No; i++)
      {
	if (expt->data[i].theta < PI/2.0) No_forward++;
	else No_backward++;
      }

      //Calculate poisson part of the log-likelihood (forward)
       PL = +Ne_tot_forward - No_forward*log(Ne_tot_forward) + logfactNo(No_forward);

      //Calculate poisson part of the log-likelihood (backward)
      PL += +Ne_tot_backward - No_backward*log(Ne_tot_backward) + logfactNo(No_backward);

      if (isinf(PL))
      {
	return 1e30;
      }

      std::cout << "No_forward = " << No_forward << ";  Ne_forward = " << Ne_tot_forward << std::endl;
      std::cout << "No_backward = " << No_backward << ";  Ne_backward = " << Ne_tot_backward << std::endl;

      //Calculate the event-by-event part
      double eventLike = 0;
      for (int i = 0; i < No; i++)
      {
	if (expt->data[i].theta < PI/2.0)
	{
	    setCurrentRate(&DMRate);
	    setCurrentVelInt(&VelInt_forwardBinned);
	    eventLike = f_S_forward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_forward;
	    setCurrentRate(&BGRate);
	    eventLike += f_BG_forward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy, &parameters)/Ne_BG_forward;
	}
	else
	{

	    setCurrentRate(&DMRate);
	    setCurrentVelInt(&VelInt_backwardBinned);
	    eventLike = f_S_backward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_backward;
	    setCurrentRate(&BGRate);
	    eventLike += f_BG_backward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy, &parameters)/Ne_BG_backward;
	}
	PL -=log(eventLike);
      }

      }
  }


    //---------------------------------------------------------------------------------------------------
    //--------Binned Momentum parametrisation: mode = 2--------------------------------------------------
    //---------------------------------------------------------------------------------------------------

    if (mode == 2)
    {
      double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);

      //Calculate poisson part of the log-likelihood
      double logfactNo = 0;

      for (int i = 0; i < No; i++)
      {
	logfactNo += log(i+1);
      }

      PL = +Ne - No*log(Ne) + logfactNo;

      //Calculate the event-by-event part
      for (int i = 0; i < No; i++)
      {
	  PL -= log((expt->exposure)*(expt->m_det)*DMRate(expt->data[i].energy,&parameters)/Ne);
      }
    }

    //---------------------------------------------------------------------------------------------------
    //--------Polynomial parametrisation: mode = 3-------------------------------------------------------
    //---------------------------------------------------------------------------------------------------
    if (mode == 3)
    {				//Use pointers to the rate functions to avoid writing all this out multiple times ----------------<<<<<<<<<<
      if (dir == 0)
      {
	setCurrentVelInt(&VelInt_isotropicPoly);

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

	     PL += PoissonLike(expt, parameters, &DMRate, No_bin,expt->bin_edges[i], expt->bin_edges[i+1] );

	   }
	  }
	else
	{

	    //Calculate expected number of events
	    double Ne = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
	    double Ne_BG = (expt->m_det)*(expt->exposure)*(N_expected(&BGRate, parameters));
	    double Ne_tot = Ne+Ne_BG;

	    //Calculate signal/BG fractions
	    double f_S = Ne/Ne_tot;
	    double f_BG = Ne_BG/Ne_tot;

	      //Calculate poisson part of the log-likelihood
	      PL = +Ne_tot - No*log(Ne_tot) + logfactNo(No);

	      //Calculate event-by-event part
	      double eventLike = 0;
	    for (int i = 0; i < No; i++)
	      {
		setCurrentRate(&DMRate);
		eventLike = f_S*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne;
		setCurrentRate(&BGRate);
		eventLike += f_BG*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_BG;
		PL -= log(eventLike);
	      }
	}

      }
      else if (dir == 1)
      {
	//Calculate expected number of events
	setCurrentVelInt(&VelInt_forwardPoly);
	double Ne_forward = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);
	setCurrentVelInt(&VelInt_backwardPoly);
	double Ne_backward = (expt->m_det)*(expt->exposure)*N_expected(&DMRate, parameters);

	//Expected background assuming isotropy
	double Ne_BG_forward = 0.5*(expt->m_det)*(expt->exposure)*N_expected(&BGRate, parameters);
	double Ne_BG_backward = Ne_BG_forward;


	double Ne_tot_forward = Ne_forward + Ne_BG_forward;
	double Ne_tot_backward = Ne_backward + Ne_BG_backward;

	//Calculate signal/BG fractions
	double f_S_forward = Ne_forward/(Ne_tot_forward);
	double f_BG_forward = 1 - f_S_forward;
	double f_S_backward = Ne_backward/(Ne_tot_backward);
	double f_BG_backward = 1 - f_S_backward;


	int No_forward = 0;
	int No_backward = 0;

	for (int i = 0; i < No; i++)
	{
	  if (expt->data[i].theta < PI/2.0) No_forward++;
	  else No_backward++;
	}

	if ((Ne_backward < 0)||(Ne_forward < 0)) return 1e30;


	//Calculate poisson part of the log-likelihood (forward)
	PL = +Ne_tot_forward - No_forward*log(Ne_tot_forward) + logfactNo(No_forward);

	//Calculate poisson part of the log-likelihood (backward)
	PL += +Ne_tot_backward - No_backward*log(Ne_tot_backward) + logfactNo(No_backward);



	//Calculate the event-by-event part
	double eventLike = 0;
	for (int i = 0; i < No; i++)
	{
	  if (expt->data[i].theta < PI/2.0)
	  {
	      setCurrentRate(&DMRate);
	      setCurrentVelInt(&VelInt_forwardPoly);
	      eventLike = f_S_forward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_forward;
	      setCurrentRate(&BGRate);
	      eventLike += f_BG_forward*(expt->exposure)*(expt->m_det)*0.5*convolvedRate(expt->data[i].energy, &parameters)/Ne_BG_forward;
	  }
	  else
	  {
	      setCurrentRate(&DMRate);
	      setCurrentVelInt(&VelInt_backwardPoly);
	      eventLike = f_S_backward*(expt->exposure)*(expt->m_det)*convolvedRate(expt->data[i].energy,&parameters)/Ne_backward;
	      setCurrentRate(&BGRate);
	      eventLike += f_BG_backward*(expt->exposure)*(expt->m_det)*0.5*convolvedRate(expt->data[i].energy, &parameters)/Ne_BG_backward;
	  }
	  PL -= log(eventLike);
	  //std::cout << PL << std::endl;
	}

      }


    }

  if (isinf(PL))
  {
    PL = 1e30;
  }
  //std::cout << PL << std::endl;
  return PL;

}

