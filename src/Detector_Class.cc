
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>

#include "DMUtils.h"
#include "Detector_Class.h"
#include "Neutrinos.h"
#include "ParamSet_Class.h"

#include "gsl/gsl_sf.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

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
	char file_FF[30];
	EFTparams.resize(N_isotopes);
	maxpow.resize(N_isotopes);
	
    for (int i = 0; i < N_isotopes; i++)
    {
      sprintf(numstr, "%d", i+1);
      m_n.push_back(read_param_double(&file, "m_"+std::string(numstr)));
      N_p.push_back(read_param_int(&file, "N_p_"+std::string(numstr)));
      N_n.push_back(read_param_int(&file, "N_n_"+std::string(numstr)));
      frac_n.push_back(read_param_double(&file, "frac_"+std::string(numstr)));
	  

	  
      //Read in nuclear response function parameters
	  sprintf(file_FF, "data/FormFactors_Z=%d_A=%d.dat", N_p[i], (N_n[i]+N_p[i]));
	  std::vector<double> temp;
	  std::ifstream ifs(file_FF);
	  double val;
	  while (ifs >> val)
	  {
		  //std::cout << val << std::endl;
		  temp.push_back(val);
	  }
	  EFTparams[i] = temp;

	  
      //Reset to start of file
      ifs.clear();
      ifs.seekg(0, std::ios::beg);
	  std::string line;
	  std::getline(ifs, line);
	  std::istringstream iss(line);
	  
	  maxpow[i] = -1;
	  while (iss >> val)
	  {
		  maxpow[i]++;
	  }

      //Read in spin parameters
      J.push_back(read_param_double(&file, "J_"+std::string(numstr)));

      if (J[i]*J[i] > 1e-6)
      {
		  /*
	Sn.push_back(read_param_double(&file, "Sn_"+std::string(numstr)));
	Sp.push_back(read_param_double(&file, "Sp_"+std::string(numstr)));

	N.push_back(read_param_double(&file, "N00_"+std::string(numstr)));
	N.push_back(read_param_double(&file, "N11_"+std::string(numstr)));

	alpha.push_back(read_param_double(&file, "alpha00_"+std::string(numstr)));
	alpha.push_back(read_param_double(&file, "alpha11_"+std::string(numstr)));

	beta.push_back(read_param_double(&file, "beta00_"+std::string(numstr)));
	beta.push_back(read_param_double(&file, "beta11_"+std::string(numstr)));
		  */
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
	  ifs.close();
    }


    m_det = read_param_double(&file, "m_det");
    exposure = read_param_double(&file, "exposure");
    E_min = read_param_double(&file, "E_min");
    E_max = read_param_double(&file, "E_max");
	//start_time = read_param_double(&file, "start_time");
    
      
    BG_level = read_param_double(&file, "BG_level");
    dE = read_param_double(&file,"dE");
    USE_BINNED_DATA = read_param_int(&file,"USE_BINNED_DATA");
    if (USE_BINNED_DATA == -1) USE_BINNED_DATA = 0;
    bin_width = read_param_double(&file, "bin_width");
    //tbin_width = read_param_double(&file, "tbin_width");
	
	USE_DIR = read_param_double(&file, "USE_DIR");
	if (USE_DIR == -1) USE_DIR = 0;
	
	DATA_FORMAT = read_param_int(&file,"DATA_FORMAT");
	if (DATA_FORMAT == -1) DATA_FORMAT = 1;
	
    file.close();
  }
  else std::cout << "Unable to open experimental parameter file:\t'" << filename << "'" << std::endl;




  //Calculate bin edges if required
  if (bin_width > 1e-3)
  {
	  //std::cout << "Warning [Detector_Class.cc]: binned data may be buggy..." << std::endl;
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
      
	/*
    tbin_edges.clear();
    double t = start_time;
    while (t <= start_time+exposure)
    {
        tbin_edges.push_back(t);
        t += tbin_width;
    }
	*/
    asimov_data.resize(N_Ebins,0.0);

    //for (int i = 0; i < (N_Ebins+1); i++)
    //{
    // std::cout << bin_edges[i] << "\t";
    //}
    //std::cout << std::endl;
	
	if (INCLUDE_NU)
	{
	//Load neutrino data
	A_nu = 1.0;
	ParamSet parameters(this,NULL, NULL);
	LoadFluxTable();
	//double total_N = 0.0;
	//std::cout << E_min_v(1e4, m_n[0]) << std::endl;
	//double p = N_expected(&NeutrinoRate,parameters, 1e4 - E_min_v(1e4, m_n[0]),10e4);
	//double r = N_expected(&NeutrinoRate,parameters, E_min_v(1e4, m_n[0]),10e4);
	//std::cout << "Here: " << p/r << std::endl;
	for (int i = 0; i < N_Ebins; i++)
	{
	    neutrino_data.push_back(m_det*exposure*N_expected(&NeutrinoRate,parameters, bin_edges[i], bin_edges[i+1]));
		//total_N += neutrino_data[i];
		//neutrino_data.push_back(0.0);
		//std::cout << neutrino_data[i] << std::endl;
    }
	ClearFluxTable();
}
	//std::cout << "Total number of neutrino scattering events is " <<  total_N << std::endl;


  }
//Load in the efficiency table
  int read_eff = 0;
  if (read_eff)
  {
	std::cout << "Loading in the efficiency table..." << std::endl;
	  //Read in nuclear response function parameters

		std::ifstream eff_fs("Efficiency-LUX.txt");
	if (eff_fs.is_open())
	{
		double Eval, effval;
		while (eff_fs >> Eval >> effval)
		{
			eff_E.push_back(Eval);
			eff_eta.push_back(effval);
		}
		eff_fs.close();
	
		//for (int i = 0; i < 50; i++)
		//{
		//		std::cout << eff_E[i] << "\t" << eff_eta[i] << std::endl;
		//}		
			
	}
	else
	{
		std::cout << " Warning [Detector_Class.cc]: Efficiency file not found..." << std::endl;
	}
	}
}

//------------------------------------------
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

//--------------------------------------------
void Detector::angular_bin_data(int N_ang_bins)
{
	std::cout << " Binning angular data..." << std::endl;
	std::vector<Event> data1;
	for (int i = 0; i < N_ang_bins; i++)
	{
		data1.clear();
		double theta1 = PI*(i)/N_ang_bins;
		double theta2 = PI*(i+1)*1.0/N_ang_bins;
		for (int j = 0; j < data.size(); j++)
		{
			if ((data[j].theta > theta1)&&(data[j].theta < theta2))
			{
				data1.push_back(data[j]);
			}	
		}
		std::cout << "\t k = " << (i+1) << ": " << data1.size() << " events" << std::endl;
		data_ang.push_back(data1);
	}
	std::cout << " Finished binning angular data..." << std::endl;

	
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
    //double R1 = 1.2*pow(m_n[i_isotope],(1.0/3.0));
    //s = 1;
    //double R1 = sqrt(x*x - 5*s*s);
   

    //Calculate Helm Form Factor squared
    double J1 = gsl_sf_bessel_j1(q*R1);

    F = 3*J1/(q*R1);
    F = pow(F,2)*exp(-pow((q*s),2));
	
  }

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
  //double b = 1.0*pow(m_n[i_isotope],1.0/6.0);

  double b = sqrt((41.467/(45.0*pow(m_n[i_isotope],-1.0/3.0) - 25.0*pow(m_n[i_isotope],-2.0/3.0))));

  double u = (q*b*q*b)/2.0;


  //----------Got rid of N[i]------

  double F = N[i]*((1-beta[i])*exp(-alpha[i]*u) + beta[i]);

  if (i_isotope==0)
  {
     //------CEFT------
     //F = exp(-u)*(0.0417889 + u*-0.111171 + pow(u,2)*0.171966 + pow(u,3)*-0.133219 + pow(u,4)*0.0633805 + pow(u,5)*-0.0178388 + pow(u,6)*0.00282476 + pow(u,7)*-2.31681e-4 + pow(u,8)*7.78223e-6 + pow(u,9)*-4.49287e-10);

     //------NijmegenII-------
     //F = exp(-2*u)*(0.0277344 + u*-0.124487 + pow(u,2)*0.328287 + pow(u,3)*-0.481399 +pow(u,4)*0.475646 + pow(u,5)*-0.285177 + pow(u,6)*0.0968193 +pow(u,7)*-0.0170957 + pow(u,8)*0.00123738);

  }
  if (i_isotope==1)
{
    //-----CEFT------
    //F = exp(-u)*(0.054731 + pow(u,1)*-0.146897 + pow(u,2)*0.182479 + pow(u,3)*-0.128112 + pow(u,4)*0.0539978 + pow(u,5)*-0.0133335 + pow(u,6)*0.00190579 + pow(u,7)*-1.48373e-4 + pow(u,8)*5.11732e-6 + pow(u,9)*-2.06597e-8);


     //----NijmegenII------
     //F = exp(-2*u)*(0.046489 +u*-0.225507 +pow(u,2)*0.499045 + pow(u,3)*-0.62243 + pow(u,4)*0.46361 + pow(u,5)*-0.20375 + pow(u,6)*0.0510851 + pow(u,7)*-0.00670516 + pow(u,8)*0.00035659);

}
  //std::cout << F << std::endl;

  return pow(F,1);
}

double Detector::EFT_formfactor(double E, int i_op, int i_isotope, int i_component)
{
    //Define conversion factor from amu-->keV
    double amu = 931.5*1000;

    //Convert recoil energy to momentum transfer q in keV
    double  q = sqrt(2*m_n[i_isotope]*amu*E);

    //Convert q into fm^-1
    q *= (1e-12/1.97e-7);
    //q *= 1e-6;

    //NB: b is in fm
    //double b = 1.0*pow(m_n[i_isotope],1.0/6.0);
    double b = sqrt((41.467/(45.0*pow(m_n[i_isotope],-1.0/3.0) - 25.0*pow(m_n[i_isotope],-2.0/3.0))));
    double y = (q*b*q*b)/4.0;
	
	if (i_op == 7)
	{	
	    double Fnn = 0.000607 - 0.00136*y + 0.000266*pow(y,2) + 0.000550*pow(y,3) + 0.0000997*pow(y,4);
	    double Fpp = 1.81 - 4.85*y + 4.88*pow(y,2) - 2.18*pow(y,3) + 0.364*pow(y,4);
	    double Fpn = -0.0331 + 0.0815*y - 0.0511*pow(y,2) - 0.00142*pow(y,3) + 0.00602*pow(y,4);
	    double F = 0.25*(Fnn + Fpp + 2*Fpn)*exp(-2*y);
		if (F < 0) std::cout << "Warning: F < 0..." << std::endl;
		return F;
		
		
	}
	else std::cout << "Parameter i_op = " << i_op << " not supported..." << std::endl;
	
	
	
}

double Detector::SI_enhancement(int i_isotope)
{	
	return  m_n[i_isotope]*m_n[i_isotope];
	//return N_p[i_isotope]*N_p[i_isotope];
}

double Detector::SD_enhancement(int i_isotope) //Watch out there is dependence on ap and an in here...check plus and minus signs...
{
  double a_p = +1;
  double a_n = +1;

  //return (8.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);


  //THIS ONE IS THE CORRECT ONE!!!
  //return (4.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);


  //return (4.0/3.0)*((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n+a_p),2);

  return (16.0*PI/3.0)*(1.0/(2*J[i_isotope]+1));

   //THIS MATCHES ARINA ET AL!!!!
   //return ((J[i_isotope]+1.0)/J[i_isotope])*pow((a_n*Sn[i_isotope]+a_p*Sp[i_isotope]),2);
}

//------------EFT response functions-----------------------------------
//------------Currently using isospin = 0 only-------------------------
double Detector::F_EFT(double E, int i_isotope, int i_component, int i_op)
{	
    //Define conversion factor from amu-->keV
    double amu = 931.5*1000;

    //Convert recoil energy to momentum transfer q in keV
    double  q = sqrt(2*m_n[i_isotope]*amu*E);

    //Convert q into fm^-1
    q *= (1e-12/1.97e-7);
    //q *= 1e-6;

    //NB: b is in fm
    //double b = 1.0*pow(m_n[i_isotope],1.0/6.0);
    double b = sqrt((41.467/(45.0*pow(m_n[i_isotope],-1.0/3.0) - 25.0*pow(m_n[i_isotope],-2.0/3.0))));
    double y = (q*b*q*b)/4.0;
	
	std::vector<double> vals = EFTparams[i_isotope];
	//int n_rows = vals.size()/(4*5);
	int n_rows = vals.size()/(4*maxpow[i_isotope]);
	
    for (int i = 0; i < vals.size(); i++)
	{
		//std::cout << vals[i] << std::endl;
		
	}
	//std::cout << std::endl;
	//std::cout << std::endl;

    //std::cout << "A\t" << m_n[i_isotope] << std::endl;
    //std::cout << "Operator\t" << i_op << std::endl;

    int ny = maxpow[i_isotope] +1;

	//std::cout << ny << std::endl;
	double coeff = 0;
	double res = 0;
	
	//std::cout << std::endl;
	for (int i = 0; i < ny; i++)
	{
		coeff = 0;
		//Need to be careful - what have I defined here!
		
		//Need to fix the stuff that's in here...it doesn't make sense...
		//Code in a_p/a_n
		
		int ind = 4*ny*i_op + ny*i_component + i;
		coeff = vals[ind];
	    //std::cout << ind << "\t" << coeff << std::endl;
		/*
		for (int j = 0; j < 4; j++)
		{
			if (j==0 || j==3)
			{
		        coeff += vals[4*ny*i_op + ny*j + i];
		    }
			if (j==1 || j==2)
			{
		        coeff += vals[4*ny*i_op + ny*j + i];
		    }
		
			//std::cout << (4*ny*i_op + ny*j + i) << "\t" << vals[4*ny*i_op + ny*j + i] << std::endl;
	    }
		*/
		res += coeff*pow(y,i);
	}
	//std::cout << res << std::endl;
	//std::cout << std::endl;
	//res = 1;
	return res*exp(-2*y);
}

double Detector::F_M(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 0);}

double Detector::F_Sigma1(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 1);}

double Detector::F_Sigma2(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 2);}

double Detector::F_Delta(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 3);}

double Detector::F_Phi2(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 4);}

double Detector::F_MPhi2(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 5);}

double Detector::F_Sigma1Delta(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 6);}

double Detector::F_SD(double E, int i_isotope, int i_component)
	{return F_EFT(E, i_isotope, i_component, 7);}

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

double Detector::efficiency(double E)
{
	if (E <= 1.1) return 0;
	if (E >= 60.0) return 0;
	
	double dlE = (log10(70) - log10(1.1))/49.0;
	
	double lE = log10(E);
	int i = floor((lE - log10(1.1))/dlE);
	
	double leta = eff_eta[i] + (lE - eff_E[i])*(eff_eta[i+1]-eff_eta[i])/dlE;
	
	//std::cout << pow(10.0,leta) << std::endl;
	//return 0.29*pow(log10(E-6.0),0.9);
	//return 1;
	return pow(10.0,leta);
}

void Detector::load_data(std::string filename)
{
  double E_value;
  double theta_value;
  double phi_value;
  
  double count = 0;
  
  double x_value, y_value, z_value;
  //Watch out for the new line!!!

  data.clear();

  std::ifstream file (filename.c_str());
  if (file.is_open())
  {  	  
    while ( (file.good())&&(!file.eof()))
    {
		if (DATA_FORMAT == 1)
		{
      	 	file >> E_value >> theta_value >> phi_value;
			data.push_back(Event(E_value,theta_value,phi_value));
  	  	}
		else if (DATA_FORMAT == 2)
		{
			file >> x_value >> y_value >> z_value;
			data.push_back(Event(0,0,0));
			int i = data.size()-1;
			data[i].E_x = x_value;
			data[i].E_y = y_value;
			data[i].E_z = z_value;
			data[i].Vector_to_Angles();
			//data[i].displayEvent();
		}
		count++;
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
