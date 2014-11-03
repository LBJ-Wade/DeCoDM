#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


#include "DMUtils.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//File scope vairables
const gsl_rng_type * T;
 gsl_rng * r;
int USE_ASIMOV;

//--------Function Prototypes-----------

void printEvents(int No, double exposure, std::string filename);
void printAsimovEvents(double Ne, std::string filename);

extern "C" { void dsinterface_nevents_( double*, double*, double*, double*, double*, double*, double* ); }
extern "C" { void dsinterface_init_( double*, int*, int*);}
//--------Main - Event Generator---------
int main(int argc, char *argv[])
{
  //Read in command line arguments
  if (argc < 4)
  {
    std::cout << "Not enough arguments - require m_x (Gev) AND sigma_SI (cm^2) AND sigma_SD (cm^2) [and optional USE_ASIMOV(0/1)]" << std::endl;
    return 0;
  }

  if (argc == 5)
  {
    USE_ASIMOV = atoi(argv[4]);
  }

  if (argc > 5)
  {
   std::cout << "Too many arguments!" << std::endl;
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


  //Read in decay channel parameters...
  int decay_channel;
  std::ifstream paramfile ("params.ini");
  if (paramfile.is_open())
    {
      decay_channel = read_param_int(&paramfile, "decay_channel");
      paramfile.close();
    }
  else std::cout << "Unable to open parameter file:\t'" << "params.ini" << "'" << std::endl;


  //Read in speed distribution parameters

  double Ne = 0;
  double result = 0;
  //Open file for reading in distribution parameters
  std::ifstream file ("dist.txt");
  if (file.is_open())
  {
    //Read in parameter values
    int N_dist = read_param_int(&file, "N_dist");
    double rho = read_param_double(&file, "rho_x");
    double v_rms;
    double v_lag[3];
    double v_lag_av;
    double fraction;
    char numstr[21]; // enough to hold all numbers up to 64-bits


    //Initialise DarkSUSY
    dsinterface_init_(&rho, &N_dist, &decay_channel);


    for (int i = 0; i < N_dist; i++)
    {
      sprintf(numstr, "%d", i+1);
      fraction = read_param_double(&file, "fraction"+std::string(numstr));
      read_param_vector(&file, "v_lag"+std::string(numstr),v_lag, 3);
      v_rms = read_param_double(&file, "v_rms" + std::string(numstr));
      v_lag_av = sqrt(v_lag[0]*v_lag[0] + v_lag[1]*v_lag[1] + v_lag[2]*v_lag[2]);

      dsinterface_nevents_(&m_x, &sigma_SI, &sigma_SD, &fraction,  &v_rms, &v_lag_av, &result);
      std::cout << result << std::endl;
      Ne += result;
    }
    file.close();
  }
  else std::cout << "Unable to open distribution parameter file:\t'" << "dist.txt" << "'" << std::endl;


  int No = gsl_ran_poisson(r, Ne);
  if (USE_ASIMOV)
  {
    No = round(Ne);
    std::cout << "NB: Generating Asimov Data..." << std::endl;
  }

  std::cout << "# events expected at IceCube-86:\t" << Ne << std::endl;
  std::cout << "# events observed at IceCube-86:\t" << No << std::endl;

  printEvents(No,77760000.0, "ICevents/myevents.dat");
  printAsimovEvents(Ne, "ICevents/myevents_Asimov.dat");
  //free(experiments);

   gsl_rng_free (r);


  return 0;
}


//--------Function Definitions-----------


//----------Use a print/display procedure from the Event thingy...

void printEvents(int No, double exposure, std::string filename)
{
  //Save events to filename

  //Output file to output events to
  std::ofstream outputfile;
  outputfile.open (filename.c_str());
  outputfile << "#This file contains the following simulated data for IceCube-86:" << std::endl;
  outputfile << "# 1. Exposure: total live time (seconds)" << std::endl;
  outputfile << "# 2. Events: individual event information" << std::endl;
  outputfile << "#" << std::endl;
  outputfile << "#The event section has the following column content:" << std::endl;
  outputfile << "#[En]         event number" << std::endl;
  outputfile << "#[E]          dummy Nchan of event (do not use this file for spectral analysis)" << std::endl;
  outputfile << "#[cos(phi)]   reconstructed cosine (angle of incoming neutrino relative to Sun)" << std::endl;
  outputfile << "#[phi uncert] paraboloid sigma in degrees" << std::endl;
  outputfile << "#" << std::endl;
  outputfile << "###--Exposure--" << std::endl;
  outputfile << "[t]\t" <<   exposure << std::endl;
  outputfile << "###--Events--" << std::endl;

  for (int i = 0; i < No; i++)
  {
    outputfile << "[En]\t" << i << std::endl;
    outputfile << "[E]  10" << std::endl;
    outputfile << "[cos(phi)]  1.0" << std::endl;
    outputfile << "[phi uncert]  2.49472" << std::endl;
    outputfile << "###" << std::endl;
  }

  //Close output file
  outputfile.close();
}

void printAsimovEvents(double Ne, std::string filename)
{
  //Save events to filename

  //Output file to output events to
  std::ofstream outputfile;
  outputfile.open (filename.c_str());
  outputfile << Ne;

  //Close output file
  outputfile.close();
}
