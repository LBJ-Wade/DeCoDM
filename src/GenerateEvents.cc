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

#include "Distributions.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//File scope vairables
const gsl_rng_type * T;
 gsl_rng * r;

//--------Function Prototypes-----------

void generateEvents(Detector* expt,double m_x, double sigma_SI, double sigma_SD);
void printEvents(std::vector<Event> data, std::string filename);

double generateBGEvents(Detector* expt);
void generateBGEvents_Asimov(Detector* expt);

void calcRotationMatrix(double* rot_matrix, double* v_lag);
void rotateEvent(double* theta, double* phi, double* rot_matrix);

//--------Main - Event Generator---------
int main(int argc, char *argv[])
{
  //Read in command line arguments
  if (argc != 4)
  {
     std::cout << "Not enough arguments - require: m_x (Gev); sigma_SI (cm^2); sigma_SD (cm^2)" << std::endl;
    return 0;
  }

  //DM mass
  double m_x = atof(argv[1]);

  //DM SI Cross-section
  double sigma_SI = atof(argv[2]);

  //DM SD Cross-section
  double sigma_SD = atof(argv[3]);

  //Initialise RNG
  gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_rng_set(r,(unsigned)time(NULL));

  //Load global parameters
  load_params("params.ini");

  std::vector<Detector> experiments;

  //------------------------------------------------------------------------------------------------
  //Check the existence of the experiment files------------------------------------------------------
  //-------------------------------------------------------------------------------------------------

  char numstr[21]; // enough to hold all numbers up to 64-bits

  for (int i = 0; i < N_expt; i++)
  {
      sprintf(numstr, "%d", i+1);
      experiments.push_back(Detector(expt_folder + "Experiment"+std::string(numstr)+".txt"));
      std::cout << "*********************************************" << std::endl;
      std::cout << "***** " << "Experiment"+std::string(numstr) << " *******************" << std::endl;
      std::cout << "*********************************************" << std::endl;
      experiments[i].displayParameters();
      generateEvents(&(experiments[i]), m_x, sigma_SI,sigma_SD);

      experiments[i].print_data(events_folder + "Events"+std::string(numstr)+".txt");
      experiments[i].print_asimov_data(events_folder + "Asimov_Events"+std::string(numstr)+".txt");
  }


  //Write a file explaining what inputs were used
  std::ofstream file("input.txt");
  if (file.is_open())
    {
      file << "m_x/GeV\t" << m_x << std::endl;
      file << "sigma_SI/cm^2\t" << sigma_SI << std::endl;
      file << "sigma_SD/cm^2\t" << sigma_SD << std::endl;

      file.close();
    }
  else std::cout << "Unable to open file input.txt" << std::endl;


  //free(experiments);

   gsl_rng_free (r);

  return 0;
}


//--------Function Definitions-----------

void generateEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD)
{
  //Generate events and store as a vector in data

  double Ne = 0;
  double Ne_BG = 0;

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
    Ne = scaling*expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));

    int No = gsl_ran_poisson(r,Ne);

    setCurrentRate(&DMRate);

    //Initialise rotation matrix
    double rot_matrix[9];
    //calcRotationMatrix(rot_matrix,astro.v_lag);

    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);

      double p_max = convolvedRate(expt->E_min,&parameters);
      double p = convolvedRate(E,&parameters);

      if (gsl_rng_uniform(r) < p/p_max)
      {
	double y = -10;
	while ((y <= -1)||(y >= 1)) //Is this the correct way of generating things...?
	{
	  //-----------------double check that this is the correct distribution-----------------------
	  //y = ((v_min(E,expt->m_n[0],m_x)/astro.v_lag) + gsl_ran_gaussian(r,astro.v_rms/astro.v_lag));
          y = gsl_ran_flat(r,-1,1);
	}
	double theta = acos(y);

	//Rotate into correct direction - FIX THIS IT's NOT ACTUALLY ROTATING!!!
	//rotateEvent(&theta,&phi,rot_matrix);

	expt->data.push_back(Event(E,theta,phi));
	N++;
      }

    }

   //Display signal event numbers
    std::cout << "Signal:\t\t # expected = " << Ne << "; # observed = " << expt->No() << std::endl;

    //Add BG events
    int No_BG = expt->No();
    Ne_BG = generateBGEvents(expt);
    No_BG = expt->No() - No_BG;

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
    }

    //Display BG event numbers
    std::cout << "Background:\t # expected = " << Ne_BG << "; # observed = " << No_BG << std::endl;

    //Display total event numbers
    std::cout << "Total:\t\t # expected = " << Ne_BG+Ne << "; # observed = " << expt->No() << std::endl;
}


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


