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

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_sort_double.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//File scope vairables
const gsl_rng_type * T;
 gsl_rng * r;

//--------Function Prototypes-----------

void generateEvents(Detector* expt,double m_x, double sigma_SI, double sigma_SD);
void printEvents(std::vector<Event> data, std::string filename);

double generateMaxwellEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc);
void generateMaxwellEvents_Asimov(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc);

double generateLisantiEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double v0, double v_esc, double k);
void generateLisantiEvents_Asimov(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double v0, double v_esc, double k);

double generateRandomEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD);

double generateBGEvents(Detector* expt);
void generateBGEvents_Asimov(Detector* expt);

void printSpectrum(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc);

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

  std::cout << "NB: Directional data are not accurate..." << std::endl << std::endl;

  //Initialise experimental parameters

  char numstr[21]; // enough to hold all numbers up to 64-bits


  //Detector* experiments;
  //experiments = (Detector*)malloc(N_expt*sizeof(Detector));

  std::vector<Detector> experiments;

  //------------------------------------------------------------------------------------------------
  //Check the existence of the experiment files------------------------------------------------------
  //-------------------------------------------------------------------------------------------------


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
      if (experiments[i].USE_BINNED_DATA) experiments[i].print_asimov_data(events_folder + "Asimov_Events"+std::string(numstr)+".txt");
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

 std::string dist_type;

 char numstr[21]; // enough to hold all numbers up to 64-bits

 double Ne = 0;
 double Ne_BG = 0;



 //Open file for reading in distribution parameters
  std::ifstream file ("dist.txt");
  if (file.is_open())
  {

    dist_type = read_param_string(&file, "dist_type");


    if (dist_type == "maxwell")
    {
      std::cout << "Using 'maxwell' type distribution..." << std::endl;

      //Distribution parameters for 'maxwell'
      int N_dist;
      double fraction;
      double sigma_v;
      double v_lag[3];
      double v_esc;

      //Read in parameter values
      N_dist = read_param_int(&file, "N_dist");
      v_esc = read_param_double(&file, "v_esc");

      for (int i = 0; i < N_dist; i++)
      {
	sprintf(numstr, "%d", i+1);
	fraction = read_param_double(&file, "fraction"+std::string(numstr));
	read_param_vector(&file, "v_lag"+std::string(numstr),v_lag);
	sigma_v = read_param_double(&file, "sigma_v" + std::string(numstr));

	//std::cout << sigma_v << '\t' << v_lag[2] << std::endl;
	Ne += generateMaxwellEvents(expt, m_x, fraction*sigma_SI, fraction*sigma_SD, v_lag, sigma_v, v_esc);
	if (expt->USE_BINNED_DATA) generateMaxwellEvents_Asimov(expt, m_x, fraction*sigma_SI, fraction*sigma_SD, v_lag, sigma_v, v_esc);

      }
    }
    else if (dist_type == "lisanti")
    {
      std::cout << "Using 'lisanti' type distribution..." << std::endl;

      //Distribution parameters for 'lisanti'
      double v0;
      double v_esc;
      double k;

      //Read in parameter values
      v0 = read_param_double(&file, "v0");
      v_esc = read_param_double(&file, "v_esc");
      k = read_param_double(&file, "k");

      Ne = generateLisantiEvents(expt, m_x, sigma_SI, sigma_SD, v0, v_esc, k);
      if (expt->USE_BINNED_DATA) generateLisantiEvents_Asimov(expt, m_x, sigma_SI, sigma_SD, v0, v_esc, k);

    }
    else if (dist_type == "random")
      {
          std::cout << "Using 'random' type distribution..." << std::endl;
	  Ne = generateRandomEvents(expt, m_x, sigma_SI, sigma_SD);

      }
    else
    {
     std::cout << "dist_type '" <<  dist_type << "' is not valid. Exiting..." << std::endl;
     exit (EXIT_FAILURE);
    }

    //Display signal event numbers
    std::cout << "Signal:\t\t # expected = " << Ne << "; # observed = " << expt->No() << std::endl;

    //Add BG events
    int No_BG = expt->No();
    if (expt->USE_BINNED_DATA) generateBGEvents_Asimov(expt);
    Ne_BG = generateBGEvents(expt);

    No_BG = expt->No() - No_BG;

    //Display BG event numbers
    std::cout << "Background:\t # expected = " << Ne_BG << "; # observed = " << No_BG << std::endl;

    //Display total event numbers
    std::cout << "Total:\t\t # expected = " << Ne_BG+Ne << "; # observed = " << expt->No() << std::endl;
    file.close();
  }
  else std::cout << "Unable to open distribution parameter file:\t'" << "dist.txt" << "'" << std::endl;

  std::cout << std::endl;



}

double generateRandomEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD)
{
    //Arrange parameters in an array
    double params[3+10];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);

    double g[11];
    for (int i = 1; i < 10; i++)
    {
     g[i] = gsl_ran_flat(r, 0, 1);
    }
    g[0] = 0;
    g[10] = 1;
    gsl_sort(g, 1, 11);
    double h[10];
    for (int i = 0; i < 10; i++)
    {
     h[i] = g[i+1] - g[i];
     params[i+3] = h[i];
    }
    N_terms = 10;

  std::ofstream file("rand.txt");
  if (file.is_open())
    {
      for (int i = 0; i < 10; i++)
      {
	file << h[i] << std::endl;
      }
      file.close();
    }
  else std::cout << "Unable to open file rand.txt" << std::endl;

    ParamSet parameters(expt,params);

    //Calculate number of expected and observed events
    setCurrentVelInt(&VelInt_isotropicBinned);
    double Ne = expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));
    int No = gsl_ran_poisson(r,Ne);

    setCurrentRate(&DMRate);

    //Initialise rotation matrix
    //double rot_matrix[9];
    // calcRotationMatrix(rot_matrix,v_lag);

    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);

      double p_max = convolvedRate(expt->E_min,&parameters);
      double p = convolvedRate(E,&parameters);

      if (gsl_rng_uniform(r) < p/p_max)
      {
	double y = 0;
	while ((y <= -1)||(y >= 1)) //Is this the correct way of generating things...?
	{
	  //-----------------double check that this is the correct distribution-----------------------
	  //y = ((v_min(E,expt->m_n[0],m_x)/v_lag_length) + gsl_ran_gaussian(r,sigma_v/v_lag_length));
	}
	double theta = acos(y);

	//Rotate into correct direction
	//rotateEvent(&theta,&phi,rot_matrix);

	expt->data.push_back(Event(E,theta,phi));
	N++;
      }

    }

    //Return number of expected events (may be required)
    return Ne;
}


double generateMaxwellEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc)
{
    //Arrange parameters in an array
    double params[8];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);
    params[3] = v_lag[0];
    params[4] = v_lag[1];
    params[5] = v_lag[2];
    params[6] = sigma_v;
    params[7] = v_esc;

    //Calculate v_lag_length
    double v_lag_length = sqrt(pow(v_lag[0],2) + pow(v_lag[1],2) + pow(v_lag[2],2));

    ParamSet parameters(expt,params);

    //Calculate number of expected and observed events
    setCurrentVelInt(&VelInt_maxwell);
    double Ne = expt->m_det*expt->exposure*(N_expected(&DMRate, parameters));
    int No = gsl_ran_poisson(r,Ne);

    setCurrentRate(&DMRate);

    //Initialise rotation matrix
    double rot_matrix[9];
    calcRotationMatrix(rot_matrix,v_lag);

    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);

      double p_max = convolvedRate(expt->E_min,&parameters);
      //double p_max = convolvedRate(15.0*reduced_m_GeV(expt->m_n[0], m_x)/(expt->m_n[0]*931.5e-3),&parameters);
      //if (p_max1 > p_max) p_max = p_max1;
      double p = convolvedRate(E,&parameters);

      if (gsl_rng_uniform(r) < p/p_max)
      {
	double y = -10;
	while ((y <= -1)||(y >= 1)) //Is this the correct way of generating things...?
	{
	  //-----------------double check that this is the correct distribution-----------------------
	  y = ((v_min(E,expt->m_n[0],m_x)/v_lag_length) + gsl_ran_gaussian(r,sigma_v/v_lag_length));
	}
	double theta = acos(y);

	//Rotate into correct direction
	rotateEvent(&theta,&phi,rot_matrix);

	expt->data.push_back(Event(E,theta,phi));
	N++;
      }

    }

    //Return number of expected events (may be required)
    return Ne;
}


void generateMaxwellEvents_Asimov(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc)
{
  //Only works with energies at the minute

  //Arrange parameters in an array
    double params[8];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);
    params[3] = v_lag[0];
    params[4] = v_lag[1];
    params[5] = v_lag[2];
    params[6] = sigma_v;
    params[7] = v_esc;


    ParamSet parameters(expt,params);

    for (int i = 0; i < expt->N_Ebins; i++)
    {
      //Calculate number of expected and observed events
      setCurrentVelInt(&VelInt_maxwell);
      double Ne = expt->m_det*expt->exposure*(N_expected(&DMRate, parameters,expt->bin_edges[i], expt->bin_edges[i+1]));
      expt->asimov_data[i] += Ne;
    }
}



double generateLisantiEvents(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double v0, double v_esc, double k)
{
    //Arrange parameters in an array
    double params[6];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);
    params[3] = v0;
    params[4] = v_esc;
    params[5] = k;


    ParamSet parameters(expt,params);

    //std::cout << "------NB: Lisanti directional data is not accurate -------------" << std::endl;



    double Norm = 0;

     double velParams[3];
     velParams[0] = params[3];
     velParams[1] = params[4];
     velParams[2] = params[5];

     Norm = Lisanti_norm(velParams);

    //Calculate number of expected and observed events
    setCurrentVelInt(&VelInt_Lisanti);
    double Ne = expt->m_det*expt->exposure*(N_expected(&DMRate, parameters))/Norm;
    int No = gsl_ran_poisson(r,Ne);


    setCurrentRate(&DMRate);


    for (int N = 0; N < No; )
    {
      double E = gsl_ran_flat(r,expt->E_min,expt->E_max);

      double phi = gsl_ran_flat(r,0,2*PI);
      double theta = acos(gsl_ran_flat(r, -1, 1));

      double p_max = convolvedRate(expt->E_min,&parameters);
      double p = convolvedRate(E,&parameters);

      if (gsl_rng_uniform(r) < p/p_max)
      {
	expt->data.push_back(Event(E,theta,phi));
	N++;
      }

    }

    //Return number of expected events (may be required)
    return Ne;
}


void generateLisantiEvents_Asimov(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double v0, double v_esc, double k)
{
  //Only works with energies at the minute

 //Arrange parameters in an array
    double params[6];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);
    params[3] = v0;
    params[4] = v_esc;
    params[5] = k;

    ParamSet parameters(expt,params);

    double Norm = 0;

     double velParams[3];
     velParams[0] = params[3];
     velParams[1] = params[4];
     velParams[2] = params[5];

     Norm = Lisanti_norm(velParams);
     //std::cout << Norm << std::endl;

     setCurrentVelInt(&VelInt_Lisanti);

    for (int i = 0; i < expt->N_Ebins; i++)
    {
      //Calculate number of expected and observed events
      double Ne = expt->m_det*expt->exposure*(N_expected(&DMRate, parameters,expt->bin_edges[i], expt->bin_edges[i+1]))/Norm;
      expt->asimov_data[i] += Ne;
    }
}


double generateBGEvents(Detector* expt)
{

  //std::cout << expt->E_min << "\t" << expt->E_max << std::endl;
  ParamSet parameters(expt,NULL);
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
  ParamSet parameters(expt,NULL);
  setCurrentRate(&BGRate);


    for (int i = 0; i < expt->N_Ebins; i++)
    {
      //Calculate number of expected and observed events
      double Ne = expt->m_det*expt->exposure*(N_expected(&BGRate, parameters, expt->bin_edges[i], expt->bin_edges[i+1]));
      expt->asimov_data[i] += Ne;
    }
}

void printSpectrum(Detector* expt, double m_x, double sigma_SI, double sigma_SD, double* v_lag, double sigma_v, double v_esc)
{
  //Arrange parameters in an array
    double params[8];

    //Theoretical parameters
    params[0] = log10(m_x);
    params[1] = log10(sigma_SI);
    params[2] = log10(sigma_SD);
    params[3] = v_lag[0];
    params[4] = v_lag[1];
    params[5] = v_lag[2];
    params[6] = sigma_v;
    params[7] = v_esc;

    //Calculate v_lag_length
    double v_lag_length = sqrt(pow(v_lag[0],2) + pow(v_lag[1],2) + pow(v_lag[2],2));

    ParamSet parameters(expt,params);

        std::cout << "------NB:Maxwell directional data is not accurate -------------" << std::endl;

    //Add in a new routine which generates the background events...
    //std::cout << 15.0*reduced_m_GeV(expt->m_n[0], m_x)/(expt->m_n[0]*931.5e-3) << std::endl;
    //Calculate number of expected and observed events
    setCurrentVelInt(&VelInt_maxwell);
    double Ne = (N_expected(&DMRate, parameters));
    //int No = gsl_ran_poisson(r,Ne);

    setCurrentRate(&DMRate);

    //Write a file explaining what inputs were used
    std::ofstream file("spectrum.txt");
    if (file.is_open())
    {
      double E = 0;
      double R = 0;
       for (int i = 1; i < 1000; i++)
	{
	  E = i*0.1;
	  R = convolvedRate(E, &parameters)/Ne;
	  file << E << "\t" << R << std::endl;
	}
      file.close();
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


