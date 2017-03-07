#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>

#include "Detector_Class.h"
#include "GenerateEvents.h"
#include "DMUtils.h"


//--------Main - Event Generator---------
int main(int argc, char *argv[])
{
  //Read in command line arguments
  if (argc < 6)
  {
     std::cout << "Not enough arguments - require: m_x (Gev); lambda_p^D (GeV^-2); lambda_n^D (GeV^-2); lambda_p^Dbar (GeV^-2); lambda_n^Dbar (GeV^-2);" << std::endl;
    return 0;
  }

  //DM mass
  double m_x = atof(argv[1]);

  double N1 = 1.0;
  //Dirac couplings
  double lambda_p_D = N1*1e-9*atof(argv[2]);
  double lambda_p_Dbar = N1*1e-9*atof(argv[3]);

  //Anti-Dirac couplings
  double lambda_n_D = N1*1e-9*atof(argv[4]);
  double lambda_n_Dbar = N1*1e-9*atof(argv[5]);

 

  Particlephysics theory;
  theory.m_x = m_x;
  theory.lambda_p_D = lambda_p_D;
  theory.lambda_n_D = lambda_n_D;
  theory.lambda_p_Dbar = lambda_p_Dbar;
  theory.lambda_n_Dbar = lambda_n_Dbar;
  theory.PrintAll();
  
  //Load global parameters
  load_params("params.ini");
  
  initialiseRNG();

  loadExperiments();

  generateAllEvents(&theory);
  
  printAllData();

  //Write a file explaining what inputs were used
  std::ofstream file("input.txt");
  if (file.is_open())
    {
      file << "m_x/GeV\t" << m_x << std::endl;
      file << "lambda_p^D/GeV^-2\t" << lambda_p_D << std::endl;
      file << "lambda_n^D/GeV^-2\t" << lambda_n_D << std::endl;
      file << "lambda_p^Dbar/GeV^-2\t" << lambda_p_Dbar << std::endl;
      file << "lambda_n^Dbar/GeV^-2\t" << lambda_n_Dbar << std::endl;

      file.close();
    }
  else std::cout << "Unable to open file input.txt" << std::endl;


  //free(experiments);

   clearRNG();

  return 0;
}