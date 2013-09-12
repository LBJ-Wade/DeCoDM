
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

#include "DMUtils.h"
#include "Detector_Class.h"

#include "gsl/gsl_sf.h"

#ifndef PI
  #define PI 3.14159265358
#endif

Detector::Detector()
{

}

Detector::Detector(std::string filename)
{
  //Open detector parameter file
  std::ifstream file (filename.c_str());
  if (file.is_open())
  {
    //Read in parameter values
    m_n = read_param_double(&file, "m_n");
    m_det = read_param_double(&file, "m_det");
    exposure = read_param_double(&file, "exposure");
    E_min = read_param_double(&file, "E_min");
    E_max = read_param_double(&file, "E_max");
    Sn = read_param_double(&file, "Sn");
    Sp = read_param_double(&file, "Sp");
    zeta = read_param_double(&file, "zeta");
    N = read_param_double(&file,"N");
    alpha = read_param_double(&file, "alpha");
    beta = read_param_double(&file, "beta");
    J = read_param_double(&file, "J");
    BG_level = read_param_double(&file, "BG_level");
    dE = read_param_double(&file,"dE");
    USE_BINNED_DATA = read_param_int(&file,"USE_BINNED_DATA");
    if (USE_BINNED_DATA == -1) USE_BINNED_DATA = 0;
    bin_width = read_param_double(&file, "bin_width");
    file.close();
  }
  else std::cout << "Unable to open experimental parameter file:\t'" << filename << "'" << std::endl;

  //Calculate bin edges if required
  if (USE_BINNED_DATA)
  {
    bin_edges.clear();
   double E = E_min;
   while (E <= E_max)
   {
        bin_edges.push_back(E);
	E += bin_width;
   }

    N_Ebins = bin_edges.size() -1;
    binned_data.clear();
    binned_data.resize(N_Ebins,0);

    asimov_data.resize(N_Ebins,0.0);

    //for (int i = 0; i < (N_Ebins+1); i++)
    //{
    // std::cout << bin_edges[i] << "\t";
    //}
    //std::cout << std::endl;

  }

}

//Count data in each bin
double Detector::bin_data()
{
   for (int i = 0; i < data.size(); i++)
   {
      int j = 0;

      while ((data[i].energy < bin_edges[j])||(data[i].energy > bin_edges[j+1]))
      {
	j++;
      }
      binned_data[j]++;
   }
}

void Detector::displayParameters()
{
 std::cout << std::endl;
 std::cout << "Detector parameters:" << std::endl;
 std::cout << "\tNuclear mass\t" << m_n << std::endl;
 std::cout << "\tDetector mass\t" << m_det << std::endl;
 std::cout << "\tExposure time\t" << exposure << std::endl;
 std::cout << "\tEnergy range\t" << E_min << "-->" << E_max << std::endl;
 std::cout << "\tEnergy resolution\t" << dE << std::endl;
 std::cout << "\tBackground level\t" << BG_level << std::endl;
 std::cout << "\tUse binned data\t" << USE_BINNED_DATA << std::endl;
 if (USE_BINNED_DATA) std::cout << "\tBin width\t" << bin_width << std::endl;
 std::cout << std::endl;

 std::cout << "Spin parameters:" << std::endl;
 std::cout << "\tProton spin <S_p>\t" << Sp << std::endl;
 std::cout << "\tNeutron spin <S_n>\t" << Sn << std::endl;
 std::cout << "\tzeta\t" << zeta << std::endl;
 std::cout << "\tN\t" << N << std::endl;
 std::cout << "\talpha\t" << alpha << std::endl;
 std::cout << "\tbeta\t" << beta << std::endl;
 std::cout << "\tJ\t" << J << std::endl;
 std::cout << std::endl;
}

void Detector::displayEvents()
{
 std::cout << No() << " Events Detected:" << std::endl;
 std::cout << "\tE (keV)\ttheta(rad)\tphi(rad)" << std::endl;
 for (int i = 0; i < No(); i++)
 {
  std::cout << "\t" << data[i].energy << "\t" << data[i].theta << "\t" << data[i].phi << std::endl;
 }
}

double Detector::SI_formfactor(double E)
{
  //Define conversion factor from amu-->keV
  double amu = 931.5*1000;

  //Convert recoil energy to momentum transfer q in keV
  double  q = sqrt(2*m_n*amu*E);

  //Convert q into fm^-1
  q *= (1e-12/1.97e-7);

  double F = 0;

  if (E < 0.001)
  {
    F = 1;
  }
  else
  {
    //Calculate Nuclear parameters in fm
    double s = 0.9;
    double a = 0.52;
    double c = (1.23*(pow(m_n,(1.0/3.0))) - 0.60);

    double R1 = sqrt(c*c + 7*PI*PI*a*a/3 - 5*s*s);

    //Calculate Helm Form Factor squared
    double J1 = gsl_sf_bessel_j1(q*R1);

    F = 3*J1/(q*R1);
    F = pow(F,2)*exp(-pow((q*s),2));
  }

//std::cout << F << std::endl;

  return F;
}

double Detector::SD_formfactor(double E)
{
  //Cerdeno, Fornasa et al (2012)

  //Define conversion factor from amu-->keV
  double amu = 931.5*1000;

  //Convert recoil energy to momentum transfer q in keV
  double  q = sqrt(2*m_n*amu*E);

  //Convert q into fm^-1
  q *= (1e-12/1.97e-7);
  //q *= 1e-6;

  //NB: b is in fm
  double b = 1.0*pow(m_n,1.0/6.0);

  double u = (q*b*q*b)/2.0;

  double F = N*((1-beta)*exp(-alpha*u) + beta);

  //Correct normalisation at zero momentum
  if (E == 0) F = 1;

  return pow(F,1);
}

double Detector::SI_enhancement()
{
 return  m_n*m_n;
}

double Detector::SD_enhancement() //Watch out there is dependence on ap and an in here...check plus and minus signs...
{
  return (16.0/3.0)*((J+1.0)/J)*pow((Sn+Sp),2);
}

//----------Input and output routines-------------------------------------------

void Detector::print_data(std::string filename)
{
  //Save events to filename

  //Output file to output events to
  std::ofstream outputfile;
  outputfile.open (filename.c_str());
  for (int i = 0; i < data.size(); i++)
  {
    data[i].printEvent(&outputfile);
  }

  //Close output file
  outputfile.close();
}

void Detector::load_data(std::string filename)
{

  double E_value;
  double theta_value;
  double phi_value;
  //Watch out for the new line!!!

  std::ifstream file (filename.c_str());
  if (file.is_open())
  {
    while ( (file.good())&&(!file.eof()) )
    {
      file >> E_value >> theta_value >> phi_value;
      data.push_back(Event(E_value,theta_value,phi_value));
    }
    file.close();
  }
  else std::cout << "Unable to open events file:\t'" << filename << "'" << std::endl;


  //Beware of this fudge-----------------------------
  data.pop_back();
}

void Detector::print_asimov_data(std::string filename)
{
  //Save asimov events to filename

  //Output file to output events to
  std::ofstream outputfile;
  outputfile.open (filename.c_str());
  for (int i = 0; i < asimov_data.size(); i++)
  {
    outputfile << bin_edges[i] << "\t" << bin_edges[i+1] << "\t" << asimov_data[i] << std::endl;
  }

  //Close output file
  outputfile.close();
}

void Detector::load_asimov_data(std::string filename)
{
  //For now bin-edges are not read in!...

  asimov_data.clear();

  double N;
  double dummy;
  //Watch out for the new line!!!

  std::ifstream file (filename.c_str());
  if (file.is_open())
  {
    while ( (file.good())&&(!file.eof()) )
    {
      file >> dummy >> dummy >> N;
      asimov_data.push_back(N);
    }
    file.close();
  }
  else std::cout << "Unable to open events file:\t'" << filename << "'" << std::endl;


  //Beware of this fudge-----------------------------
  asimov_data.pop_back();
}
