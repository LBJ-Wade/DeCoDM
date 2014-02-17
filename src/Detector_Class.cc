
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
    N_isotopes = read_param_double(&file, "N_isotopes");

    char numstr[21]; // enough to hold all numbers up to 64-bits
    for (int i = 0; i < N_isotopes; i++)
    {
      sprintf(numstr, "%d", i+1);
      m_n.push_back(read_param_double(&file, "m_"+std::string(numstr)));
      frac_n.push_back(read_param_double(&file, "frac_"+std::string(numstr)));

      //Read in spin parameters
      J.push_back(read_param_double(&file, "J_"+std::string(numstr)));

      if (J[i]*J[i] > 1e-6)
      {
	Sn.push_back(read_param_double(&file, "Sn_"+std::string(numstr)));
	Sp.push_back(read_param_double(&file, "Sp_"+std::string(numstr)));

	N.push_back(read_param_double(&file, "N00_"+std::string(numstr)));
	N.push_back(read_param_double(&file, "N11_"+std::string(numstr)));

	alpha.push_back(read_param_double(&file, "alpha00_"+std::string(numstr)));
	alpha.push_back(read_param_double(&file, "alpha11_"+std::string(numstr)));

	beta.push_back(read_param_double(&file, "beta00_"+std::string(numstr)));
	beta.push_back(read_param_double(&file, "beta11_"+std::string(numstr)));
      }
      else
      {
	//Still need to pad the arrays to make sure the indexing works
	Sn.push_back(0);
	Sp.push_back(0);
	N.push_back(0);
	N.push_back(0);
	alpha.push_back(0);
	alpha.push_back(0);
	beta.push_back(0);
	beta.push_back(0);
      }
    }


    m_det = read_param_double(&file, "m_det");
    exposure = read_param_double(&file, "exposure");
    E_min = read_param_double(&file, "E_min");
    E_max = read_param_double(&file, "E_max");

    BG_level = read_param_double(&file, "BG_level");
    dE = read_param_double(&file,"dE");
    USE_BINNED_DATA = read_param_int(&file,"USE_BINNED_DATA");
    if (USE_BINNED_DATA == -1) USE_BINNED_DATA = 0;
    bin_width = read_param_double(&file, "bin_width");
    file.close();
  }
  else std::cout << "Unable to open experimental parameter file:\t'" << filename << "'" << std::endl;

  //Calculate bin edges if required
  if (bin_width > 1e-3)
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
void Detector::bin_data()
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
 std::cout << "\tDetector mass\t" << m_det << std::endl;
 std::cout << "\tExposure time\t" << exposure << std::endl;
 std::cout << "\tEnergy range\t" << E_min << "-->" << E_max << std::endl;
 std::cout << "\tEnergy resolution\t" << dE << std::endl;
 std::cout << "\tBackground level\t" << BG_level << std::endl;
 std::cout << "\tUse binned data\t" << USE_BINNED_DATA << std::endl;
 if (USE_BINNED_DATA) std::cout << "\tBin width\t" << bin_width << std::endl;
 std::cout << std::endl;

 std::cout << "Detector composition:" << std::endl;
 for (int i = 0; i < N_isotopes; i++)
 {
    std::cout << "\t" << frac_n[i]*100 << "% - nuclear mass " << m_n[i] << " amu" << std::endl;
 }
 std::cout << std::endl;

 std::cout << "Spin parameters:" << std::endl;
 /*
 for (int i = 0; i < N_isotopes; i++)
 {
   std::cout << "------------Isotope " << (i+1) << ":" << std::endl;
    std::cout << "\tJ\t" << J[i] << std::endl;
    if (J[i] > 1e-3)
    {
      std::cout << "\tProton spin <S_p>\t" << Sp[i] << std::endl;
      std::cout << "\tNeutron spin <S_n>\t" << Sn[i] << std::endl;
      std::cout << "\tN_00\t" << N[2*i] << "\tN_11\t" << N[2*i+1] << std::endl;
      std::cout << "\talpha_00\t" << alpha[2*i] << "\talpha_11\t" << alpha[2*i+1] << std::endl;
      std::cout << "\tbeta_00\t" << beta[2*i] << "\tbeta_11\t" << beta[2*i+1] << std::endl;
    }

 }
 */
 //std::cout << "\tJ\t" << J << std::endl;
 //std::cout << "\tProton spin <S_p>\t" << Sp << std::endl;
 //std::cout << "\tNeutron spin <S_n>\t" << Sn << std::endl;
 std::cout << "\t...spin parameters not shown..." << std::endl;
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

double Detector::SI_formfactor(double E, int i_isotope)
{
  //Define conversion factor from amu-->keV
  double amu = 931.5*1000;

  //Convert recoil energy to momentum transfer q in keV
  double  q = sqrt(2*m_n[i_isotope]*amu*E);

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
    double c = (1.23*(pow(m_n[i_isotope],(1.0/3.0))) - 0.60);

    double R1 = sqrt(c*c + 7*PI*PI*a*a/3 - 5*s*s);

   

    //Calculate Helm Form Factor squared
    double J1 = gsl_sf_bessel_j1(q*R1);

    F = 3*J1/(q*R1);
    F = pow(F,2)*exp(-pow((q*s),2));
  }

//std::cout << F << std::endl;

  return F;
}

double Detector::SD_formfactor(double E, int i_isotope, int i_component)
{
  //Cerdeno, Fornasa et al (2012)

  int i = 2*i_isotope + i_component;

  //Define conversion factor from amu-->keV
  double amu = 931.5*1000;

  //Convert recoil energy to momentum transfer q in keV
  double  q = sqrt(2*m_n[i_isotope]*amu*E);

  //Convert q into fm^-1
  q *= (1e-12/1.97e-7);
  //q *= 1e-6;

  //NB: b is in fm
  double b = 1.0*pow(m_n[i_isotope],1.0/6.0);

  double u = (q*b*q*b)/2.0;


  //----------Got rid of N[i]------

  double F = N[i]*((1-beta[i])*exp(-alpha[i]*u) + beta[i]);

  if (i_isotope==0)
  {
     //------CEFT------
     F = exp(-u)*(0.0417889 + u*-0.111171 + pow(u,2)*0.171966 + pow(u,3)*-0.133219 + pow(u,4)*0.0633805 + pow(u,5)*-0.0178388 + pow(u,6)*0.00282476 + pow(u,7)*-2.31681e-4 + pow(u,8)*7.78223e-6 + pow(u,9)*-4.49287e-10);

     //------NijmegenII-------
     //F = exp(-2*u)*(0.0277344 + u*-0.124487 + pow(u,2)*0.328287 + pow(u,3)*-0.481399 +pow(u,4)*0.475646 + pow(u,5)*-0.285177 + pow(u,6)*0.0968193 +pow(u,7)*-0.0170957 + pow(u,8)*0.00123738);

  }
  if (i_isotope==1)
{
    //-----CEFT------
    F = exp(-u)*(0.054731 + pow(u,1)*-0.146897 + pow(u,2)*0.182479 + pow(u,3)*-0.128112 + pow(u,4)*0.0539978 + pow(u,5)*-0.0133335 + pow(u,6)*0.00190579 + pow(u,7)*-1.48373e-4 + pow(u,8)*5.11732e-6 + pow(u,9)*-2.06597e-8);


     //----NijmegenII------
     //F = exp(-2*u)*(0.046489 +u*-0.225507 +pow(u,2)*0.499045 + pow(u,3)*-0.62243 + pow(u,4)*0.46361 + pow(u,5)*-0.20375 + pow(u,6)*0.0510851 + pow(u,7)*-0.00670516 + pow(u,8)*0.00035659);

}
  //std::cout << F << std::endl;

  return pow(F,1);
}

double Detector::SI_enhancement(int i_isotope)
{
 return  m_n[i_isotope]*m_n[i_isotope];
}

double Detector::SD_enhancement(int i_isotope) //Watch out there is dependence on ap and an in here...check plus and minus signs...
{
  double a_p = +1;
  double a_n = +1;

  //return (8.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);


  //THIS ONE IS THE CORRECT ONE!!!
  //return (4.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);


  //return (4.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n+a_p),2);

  return (4.0*PI)*(1.0/(2*J[i_isotope]+1));

   //THIS MATCHES ARINA ET AL!!!!
   //return ((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);
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
