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
#include "DMUtils.h"
#include "EventRates.h"
#include "Astrophysics_Class.h"
#include "Particlephysics_Class.h"
#include "Neutrinos.h"
#include "Distributions.h"
#include "GenerateEvents.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_sort_double.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//File scope vairables
static const gsl_rng_type * T;
static gsl_rng * r;
 

 void loadExperiments()
 {
 	static int count = 0;
 	//Check to make sure the detector data has been loaded from file...
 	if (count == 0)
 	{
         //Load global parameters
         //load_params("params.ini");
	
	     char numstr[21]; // enough to hold all numbers up to 64-bits

	     for (int i = 0; i < N_expt; i++)
	     {
	         sprintf(numstr, "%d", i+1);
	         experiments.push_back(Detector(expt_folder + "Experiment"+std::string(numstr)+".txt"));
			 if (LOUD_MOUTH)
			 {
				 /*
		         std::cout << "*********************************************" << std::endl;
		         std::cout << "***** " << "Experiment"+std::string(numstr) << " *******************" << std::endl;
		         std::cout << "*********************************************" << std::endl;
		         experiments[i].displayParameters();
				 */
		     }
	     }
		
 	    count++;
 	}
	
 }
 
 void initialiseRNG()
 {
     //Initialise RNG
     gsl_rng_env_setup();
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);
     gsl_rng_set(r,(unsigned)time(NULL));
 }
 
 void clearRNG()
 {
 	
    gsl_rng_free (r);
	
 }



 void generateAllEvents(Particlephysics* ptheory)
 {
     for (int i = 0; i < N_expt; i++)
     {
		 if (LOUD_MOUTH)
		 {
			  char numstr[21];
			  sprintf(numstr, "%d", i+1);
			  std::cout << std::endl;
	         std::cout << "*********************************************" << std::endl;
	         std::cout << "***** " << "Experiment"+std::string(numstr) << " *******************" << std::endl;
	         std::cout << "*********************************************" << std::endl;
	         experiments[i].displayParameters();	
		 }
         generateEvents(&(experiments[i]), ptheory);
	 }
 }
 
 void generateAllEvents_Ne(double m_x, int No_target, int op)
 {
     for (int i = 0; i < N_expt; i++)
     {
         generateEvents_Ne(&(experiments[i]), m_x, No_target, op);
	 }
 }
 
 void printAllData()
 {
	 char numstr[21];
     for (int i = 0; i < N_expt; i++)
     {
		 sprintf(numstr, "%d", i+1);
         experiments[i].print_data(events_folder + "Events"+std::string(numstr)+".txt");
         if (experiments[i].USE_BINNED_DATA) experiments[i].print_asimov_data(events_folder + "Asimov_Events"+std::string(numstr)+".txt");
	 }
 }



//--------Function Definitions-----------

void generateEvents(Detector* expt, Particlephysics* ptheory)
{
	//std::cout << "This code can be rewritten for generic rates (i.e. a 'generate' dataset code and then you just set which eventrate you want to use...) - 07/11/2014" << std::endl;
	
  //Generate events and store as a vector in data
  double Ne = 0;
  double Ne_BG = 0;
  double Ne_nu = 0;

  //Empty previous events
  expt->data.clear();

  //Initialise and load in astrophysics

  Astrophysics astro;

  astro.load_params();
  //std::cout << "Astro:\t" << astro.v_lag[0] << "\t" << astro.v_rms[0] << "\t" << astro.rho_x << std::endl;
  //std::cout << "Astro:\t" << astro.dist_type << std::endl;
  //Particlephysics theory;
  //theory.m_x = m_x;
  //theory.sigma_SI = sigma_SI;
  //theory.sigma_SD = sigma_SD;

  ParamSet parameters(expt,ptheory, &astro);
 
  //Need to divide through by sum over frac_n[i] A^2
  
  //There's a slight change in here - need to include 
  //factors of mu_{\chi N}^2 - they don't *quite* cancel...
  
  //Calculating the effective DM-proton coupling
  double sigp = 0;
  double nA = 0;
  for (int i = 0; i < expt->N_isotopes; i++)
  {
	  nA += expt->frac_n[i]*pow(expt->N_p[i] + expt->N_n[i],2.0);
	  sigp += expt->frac_n[i]*pow(reduced_m_GeV(1.0, ptheory->m_x),2.0)*(pow(ptheory->lambda_p_D*expt->N_p[i] + ptheory->lambda_n_D*expt->N_n[i],2.0)
			  + pow(ptheory->lambda_p_Dbar*expt->N_p[i] + ptheory->lambda_n_Dbar*expt->N_n[i],2.0));
  }
  sigp *= (2.0/PI)*(1.973e-14*1.973e-14)/nA;
  
  //std::cout << sigp << std::endl;
  //std::cout << ptheory->lambda_p_D << std::endl;
  //std::ofstream outfile;

   //outfile.open("norms.txt", std::ios_base::app);
   //outfile << sigp << std::endl; 
  //  outfile.close();
 
  //BJK!!!
  //This should only be used for m_x = 50 or m_x = 1000
  double xsec_target
  if (ptheory->m_x < 100)
  {
	  xsec_target = 1e-46;
  }
  else if (ptheory->m_x > 100)
  {
	  xsec_target = 1e-45;
  }
  
  if (expt->m_n[0] > 120)
  //if ((expt->m_n[0] < 45)&&(expt->m_n[0] > 38))
  {
	  ptheory->lambda_p_D *= sqrt(xsec_target/sigp);
	  ptheory->lambda_n_D *= sqrt(xsec_target/sigp);
	  ptheory->lambda_p_Dbar *= sqrt(xsec_target/sigp);
	  ptheory->lambda_n_Dbar *= sqrt(xsec_target/sigp);
	  
	  	sigp = 0;
	    for (int i = 0; i < expt->N_isotopes; i++)
	    {
	  	  sigp += expt->frac_n[i]*pow(reduced_m_GeV(1.0, ptheory->m_x),2.0)*(pow(ptheory->lambda_p_D*expt->N_p[i] + ptheory->lambda_n_D*expt->N_n[i],2.0)
	  			  + pow(ptheory->lambda_p_Dbar*expt->N_p[i] + ptheory->lambda_n_Dbar*expt->N_n[i],2.0));
	    }
	    sigp *= (2.0/PI)*(1.973e-14*1.973e-14)/nA;
	  
  }
  
  //sigp = 1e-46;
  
  std::cout << "   Effective DM-proton cross section: " << sigp << " cm^2..."<<std::endl;
  
  double scaling = 1;
  if (astro.dist_type == "lisanti")
    {
       scaling = 1.0/(Lisanti_norm(&astro));
    }

    //Generate ordinary events
    Ne = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));
	
	/*
	std::ofstream outfile;

	outfile.open("Ne.txt", std::ios_base::app);
	outfile << Ne << "\t"; 
	if (expt->m_n[0] < 30) outfile << std::endl;
	outfile.close();
	*/
	
    int No = gsl_ran_poisson(r,Ne);
	//std::cout << "Warning! Rounding to nearest value!..."<<std::endl;
	//int No = round(Ne);
    
	setCurrentRate(&DMRate);

    //Initialise rotation matrix
    double rot_matrix[9];
    //calcRotationMatrix(rot_matrix,astro.v_lag);
	//std::cout << v_min(expt->E_min,expt->m_n[0],m_x) << std::endl;
    double p_max = 10.0*DMRateDirectional(0.0, PI/2.0, 0.0, &parameters);
	
    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);

	  double theta = acos(gsl_ran_flat(r,-1,1));

	  //double p = convolvedRate(E,&parameters);
	  double p = DMRateDirectional(E, theta, phi, &parameters);
	  //*-----------------Note that I'm 
	  //std::cout << p << "\t" << p_max << "\t" << p/p_max << std::endl;

      if (gsl_rng_uniform(r) < p/p_max)
      {
		expt->data.push_back(Event(E,theta,phi));
		N++;
      }

    }

	//Load neutrino flux table
	LoadFluxTable();

    int Ne_sig = expt->No();
    //Add BG events
    int No_BG = expt->No();
    Ne_BG = generateBGEvents(expt);
	No_BG = expt->No() - No_BG;
	int No_nu = expt->No();
	if (INCLUDE_NU) Ne_nu = generateNeutrinoEvents(expt);
    No_nu = expt->No() - No_nu;


    //Calculate ASIMOV data if the bin width is defined
    if ((expt->bin_width > 1e-3))
    {
	    double Ne_bin = 0;

	    for (int i = 0; i < expt->N_Ebins; i++)
		{
		  Ne_bin = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters,expt->bin_edges[i], expt->bin_edges[i+1]));
		  expt->asimov_data[i] += Ne_bin;
		}
	
		generateBGEvents_Asimov(expt);
		if (INCLUDE_NU) generateNeutrinoEvents_Asimov(expt);
    }
	  
	//Clear neutrino flux table from memory
	ClearFluxTable();

    if (LOUD_MOUTH)
	{
	    std::cout << std::endl;
	   //Display signal event numbers
	    std::cout << "Signal:\t\t # expected = " << Ne << "; # observed = " << Ne_sig << std::endl;

	    //Display BG event numbers
	    std::cout << "Background:\t # expected = " << Ne_BG << "; # observed = " << No_BG << std::endl;

	    //Display Neutrino event numbers
		if (INCLUDE_NU)
		{
		    std::cout << "CNS:\t\t # expected = " << Ne_nu << "; # observed = " << No_nu << std::endl;
		}
	    //Display total event numbers
	    std::cout << "Total:\t\t # expected = " << Ne_nu+Ne_BG+Ne << "; # observed = " << expt->No() << std::endl;
    }
}

double calcNe(Detector* expt, double m_x, double sigma_SI, double sigma_SD)
{	
  //Calculate expected number of events
	if (sigma_SI > 1e-30) return 1e10;
  //Initialise and load in astrophysics
  Astrophysics astro;
  astro.load_params();
  //std::cout << "Astro:\t" << astro.v_lag[0] << "\t" << astro.v_rms[0] << "\t" << astro.rho_x << std::endl;
  //std::cout << "Astro:\t" << astro.dist_type << std::endl;
  Particlephysics theory;
  theory.m_x = m_x;
  theory.sigma_SI = sigma_SI;
  theory.sigma_SD = sigma_SD;

  ParamSet parameters(expt,&theory, &astro);
  

  
  double scaling = 1;
  if (astro.dist_type == "lisanti")
    {
       scaling = 1.0/(Lisanti_norm(&astro));
    }

    //Generate ordinary events
    double Ne_S = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));

	double Ne_nu;
	double Ne_BG;
	//Load neutrino flux table
	if (INCLUDE_NU)
	{
	    LoadFluxTable();
		Ne_nu = generateNeutrinoEvents(expt);
		//Clear neutrino flux table from memory
		ClearFluxTable();
	}
	Ne_BG = generateBGEvents(expt);

	return Ne_S + Ne_BG + Ne_nu;
}

double calcNe_NR(Detector* expt, double m_x, int op1, int op2, int N1, int N2)
{	

  //Initialise and load in astrophysics
  Astrophysics astro;
  astro.load_params();

  Particlephysics theory;
  theory.m_x = m_x;
  
  theory.op1 = op1;
  theory.op2 = op2;
  
  theory.N1 = N1;
  theory.N2 = N2;

  ParamSet parameters(expt,&theory, &astro);
  
  double scaling = 1;
  if (astro.dist_type == "lisanti")
    {
       scaling = 1.0/(Lisanti_norm(&astro));
    }



   //Generate ordinary events
   double Ne_S = scaling*expt->m_det*expt->exposure*(N_expected(&DMRateNR, parameters));

	return Ne_S;
}


void generateEvents_Ne(Detector* expt, double m_x, int No_target, int op)
{
	//std::cout << "This code can be rewritten for generic rates (i.e. a 'generate' dataset code and then you just set which eventrate you want to use...) - 07/11/2014" << std::endl;
	
  //Generate events and store as a vector in data
  double Ne = 0;
  double Ne_BG = 0;
  double Ne_nu = 0;

  //Empty previous events
  expt->data.clear();

  //Initialise and load in astrophysics

  Astrophysics astro;


  astro.load_params();
  //std::cout << "Astro:\t" << astro.v_lag[0] << "\t" << astro.v_rms[0] << "\t" << astro.rho_x << std::endl;
  //std::cout << "Astro:\t" << astro.dist_type << std::endl;
  Particlephysics theory;
  theory.m_x = m_x;
  theory.sigma_SI = 0.0;
  theory.sigma_SD = 0.0;
  theory.sigma_O5 = 0.0;
  theory.sigma_O7 = 0.0;
  theory.sigma_O15 = 0.0;
  
  if (op == 1) theory.sigma_SI = 1e-45;
  if (op == 4) theory.sigma_SD = 1e-45;
  if (op == 5) theory.sigma_O5 = 1e-45;
  if (op == 7) theory.sigma_O7 = 1e-45;
  if (op == 15) theory.sigma_O15 = 1e-45;

  ParamSet parameters(expt,&theory, &astro);
  
  double scaling = 1;
  if (astro.dist_type == "lisanti")
    {
       scaling = 1.0/(Lisanti_norm(&astro));
    }

    //Generate ordinary events
    Ne = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));
	theory.sigma_SI *= No_target/Ne;
	theory.sigma_SD *= No_target/Ne;
	theory.sigma_O5 *= No_target/Ne;
	theory.sigma_O7 *= No_target/Ne;
	theory.sigma_O15 *= No_target/Ne;

    int No = gsl_ran_poisson(r,Ne);
	No = No_target;
    //std::cout << Ne << "\t" << No << std::endl;
    setCurrentRate(&DMRate);

    //Initialise rotation matrix
    double rot_matrix[9];
    //calcRotationMatrix(rot_matrix,astro.v_lag);
	//std::cout << v_min(expt->E_min,expt->m_n[0],m_x) << std::endl;
    double p_max = 10.0*DMRateDirectional(0.0, PI/2.0, 0.0, &parameters);
	
    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);

	  double theta = acos(gsl_ran_flat(r,-1,1));

	  //double p = convolvedRate(E,&parameters);
	  double p = DMRateDirectional(E, theta, phi, &parameters);
	  //*-----------------Note that I'm 
	  //std::cout << p << "\t" << p_max << "\t" << p/p_max << std::endl;

      if (gsl_rng_uniform(r) < p/p_max)
      {
		expt->data.push_back(Event(E,theta,phi));
		N++;
      }

    }

	//Load neutrino flux table
	LoadFluxTable();

    int Ne_sig = expt->No();
    //Add BG events
    int No_BG = expt->No();
    Ne_BG = generateBGEvents(expt);
	No_BG = expt->No() - No_BG;
	int No_nu = expt->No();
	if (INCLUDE_NU) Ne_nu = generateNeutrinoEvents(expt);
    No_nu = expt->No() - No_nu;


    //Calculate ASIMOV data if the bin width is defined
    if ((expt->bin_width > 1e-3))
    {
	    double Ne_bin = 0;

	    for (int i = 0; i < expt->N_Ebins; i++)
		{
		  Ne_bin = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters,expt->bin_edges[i], expt->bin_edges[i+1]));
		  expt->asimov_data[i] += Ne_bin;
		}
	
		generateBGEvents_Asimov(expt);
		if (INCLUDE_NU) generateNeutrinoEvents_Asimov(expt);
    }
	  
	//Clear neutrino flux table from memory
	ClearFluxTable();

    if (LOUD_MOUTH)
	{
	    std::cout << std::endl;
	   //Display signal event numbers
	    std::cout << "Signal:\t\t # expected = " << Ne << "; # observed = " << Ne_sig << std::endl;

	    //Display BG event numbers
	    std::cout << "Background:\t # expected = " << Ne_BG << "; # observed = " << No_BG << std::endl;

	    //Display Neutrino event numbers
		std::cout << "CNS:\t\t # expected = " << Ne_nu << "; # observed = " << No_nu << std::endl;

	    //Display total event numbers
	    std::cout << "Total:\t\t # expected = " << Ne_nu+Ne_BG+Ne << "; # observed = " << expt->No() << std::endl;
    }
}

//---------Generate BG events----------------------------------------------
double generateBGEvents(Detector* expt)
{

  //std::cout << expt->E_min << "\t" << expt->E_max << std::endl;
  ParamSet parameters(expt,NULL, NULL);
  double Ne_BG = expt->m_det*expt->exposure*(N_expected(&BGRate, parameters));

  int No_BG = gsl_ran_poisson(r,Ne_BG);

  double E, phi, y, theta;
  double p;

  setCurrentRate(&BGRate);
  double p_max = convolvedRate(expt->E_min,&parameters);

  for (int N = 0; N < No_BG;)
  {
     //Generate a candidate event
     E = gsl_ran_flat(r,expt->E_min,expt->E_max);
     p = convolvedRate(E,&parameters);

     //Test to see if event should be added
      if (gsl_rng_uniform(r) < p/p_max)
      {
	//Generate isotropic angular distribution
	phi = gsl_ran_flat(r,0,2*PI);
	y = gsl_ran_flat(r,-1,1);
	theta = acos(y);

	//Add event
	expt->data.push_back(Event(E,theta,phi));
	N++;
      }
   }

   //Return expected number of BG events (if needed)
   return Ne_BG;
}

void generateBGEvents_Asimov(Detector* expt)
{
  ParamSet parameters(expt,NULL, NULL);
  setCurrentRate(&BGRate);


    for (int i = 0; i < expt->N_Ebins; i++)
    {
      //Calculate number of expected and observed events
      double Ne = expt->m_det*expt->exposure*(N_expected(&BGRate, parameters, expt->bin_edges[i], expt->bin_edges[i+1]));
      expt->asimov_data[i] += Ne;
    }
}

//----------------------------------------------------------------------------------
//------------Generate neutrino background events-----------------------------------
//NB: I can potentially make this quite generic...----------------------------------
double generateNeutrinoEvents(Detector* expt)
{
	
  ParamSet parameters(expt,NULL, NULL);
  double Ne_nu = expt->m_det*expt->exposure*(N_expected(&NeutrinoRate, parameters));

  int No_nu = gsl_ran_poisson(r,Ne_nu);
  //std::cout << "Currently setting No_nu = 0" << std::endl;
  //int No_nu = 0;
  
  double E, phi, y, theta;
  double p;

  setCurrentRate(&NeutrinoRate);
  //Just so happens that the rate is constantly decreasing, so we can actually use this, reasonably
  double p_max = 1.0*convolvedRate(expt->E_min,&parameters);

  for (int N = 0; N < No_nu;)
  {
     //Generate a candidate event
     E = gsl_ran_flat(r,expt->E_min,expt->E_max);
     p = convolvedRate(E,&parameters);

     //Test to see if event should be added
      if (gsl_rng_uniform(r) < p/p_max)
      {
		//Generate isotropic angular distribution
		phi = gsl_ran_flat(r,0,2*PI);
		y = gsl_ran_flat(r,-1,1);
		theta = acos(y);

		//Add event
		expt->data.push_back(Event(E,theta,phi));
		N++;
      }
   }
   
   //Return expected number of neutrino events (if needed)
   return Ne_nu;
}

void generateNeutrinoEvents_Asimov(Detector* expt)
{
  ParamSet parameters(expt,NULL, NULL);
  setCurrentRate(&NeutrinoRate);


    for (int i = 0; i < expt->N_Ebins; i++)
    {
      //Calculate number of expected and observed events
      double Ne = expt->m_det*expt->exposure*(N_expected(&NeutrinoRate, parameters, expt->bin_edges[i], expt->bin_edges[i+1]));
      expt->asimov_data[i] += Ne;
    }
}


//------------Geometric operations to get all the angles in an x,y,z basis----------

void calcRotationMatrix(double* rot_matrix, double* v_lag)
{
  //Get length of v_lag
 double v_lag_length = sqrt(pow(v_lag[0],2) + pow(v_lag[1],2) + pow(v_lag[2],2));

 //Calculate angles for rotation
 double theta = acos(v_lag[2]/v_lag_length);
 double phi;
 if ((fabs(v_lag[0]) < 1e-20)&&(fabs(v_lag[1]) < 1e-20))
 {
   phi = 0;
 }
 else
 {
  phi = PI + atan2(v_lag[0],v_lag[1]);
 }
 //Calculate matrix elements
 rot_matrix[0] = cos(phi)*cos(theta);
 rot_matrix[1] = sin(phi);
 rot_matrix[2] = -sin(theta)*cos(phi);

 rot_matrix[3+0] = -sin(phi)*cos(theta);
 rot_matrix[3+1] = cos(phi);
 rot_matrix[3+2] = sin(phi)*sin(theta);

 rot_matrix[6+0] = sin(theta);
 rot_matrix[6+1] = 0;
 rot_matrix[6+2] = cos(theta);

}

//Does it make sense to calculate in terms of angles or keep everything in vectorised format???


void rotateEvent(double* theta, double* phi, double* rot_matrix)
{
 //Performs in place rotation of event angles for initial angles generated using v_lag pointing along z

 //Generate unit vector along initial direction
 double q_in[3];
 q_in[0] = sin(*theta)*cos(*phi);
 q_in[1] = sin(*theta)*sin(*phi);
 q_in[2] = cos(*theta);

 //Initialise unit vector along final direction
 double q_out[3] = {0,0,0};

 //Multiply by rotation matrix
 for (int i = 0; i < 3; i++)
 {
   for (int j = 0; j < 3; j++)
   {
     q_out[i] += rot_matrix[i*3 + j]*q_in[j];
   }
 }

 //Calculate final angles
  *theta = acos(q_out[2]);
  if ((fabs(q_out[0]) < 1e-20)&&(fabs(q_out[1]) < 1e-20))
  {
    *phi = 0;
  }
  else
  {
    *phi = PI + atan2(q_out[0],q_out[1]);
  }

}



//----------Use a print/display procedure from the Event thingy...


