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

  Particlephysics theory;
  theory.m_x = m_x;
  theory.sigma_SI = sigma_SI;
  theory.sigma_SD = sigma_SD;
  theory.PrintAll();
  
  //Load global parameters
  load_params("params.ini");
  
  initialiseRNG();

  loadExperiments();

  generateAllEvents(m_x, sigma_SI, sigma_SD);
  
  printAllData();

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

   clearRNG();

  return 0;
}