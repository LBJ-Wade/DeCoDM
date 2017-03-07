#ifndef DETECTORCLASS_H
#define DETECTORCLASS_H

#include <vector>
#include <string>
#include "Event_Class.h"

class Detector
{
  public:

    double N_isotopes; //Number of distinct elements/isotopes
    std::vector<double> frac_n; //Nuclear mass fraction of each element/isotope
    std::vector<double> m_n; //Nuclear mass in amu of each element/isotope

    std::vector<int> N_p; //Number of protons
    std::vector<int> N_n; //Number of neutrons
    
	std::vector<double> eff_E;
	std::vector<double> eff_eta;
	
	std::vector< std::vector<double> > EFTparams; //Parameters for EFT form factors
	std::vector<int> maxpow; //Max power of y in Form factor expansion
	
    //Spin parameters - NB: possibly introduce a nuclei class (could keep separate 'nuclear' data files and load them up when needed...)
    std::vector<double> J;
    std::vector<double> Sp;
    std::vector<double> Sn;

    std::vector<double> N;
    std::vector<double> alpha;
    std::vector<double> beta;

    double m_det; //Total detector mass in kg
    double exposure; //Effective exposure (t_exp x eff)
    double E_min;
    double E_max;
    double start_time; //Start time in days!

    //Background event levels (in events /keV/kg/day
    double BG_level;

    //Bin data and use binned likelihood analysis
    int USE_BINNED_DATA;
	
	//Format of the data ((1) - energy-angle space; or (2) - vector space)
    int DATA_FORMAT;

    //Bin width and number of bins for binned analysis
    double bin_width;
    double tbin_width;
    int N_Ebins;
    int N_tbins;
	
	//Use directional information
	int USE_DIR;

    //Energy resultion (keV)
    double dE;


    std::vector<Event> data; //List of events in detector
    std::vector<double> bin_edges; //List of energy bin edges
    std::vector<double> tbin_edges; //List of time bin edges
    std::vector<int> binned_data; //List of binned data
    std::vector<double> asimov_data; //List of asimov data in each energy/time bin
    //std::vector<double> timed_asimov_data; //List of asimov data in each time/energy bin
	double A_nu;
	std::vector<double> neutrino_data; //List of neutrino events in each bins
    double No(){return data.size();}

    void angular_bin_data(int N_ang_bins);
    std::vector<std::vector<Event> > data_ang;

   //Constructor
    Detector(); //Default constructor
    Detector(std::string); //Initialise from an experimental parameter file

    void displayParameters();
    void displayEvents();

    //Nuclear form factors
    double SI_formfactor(double E, int i_isotope);
    double SD_formfactor(double E, int i_isotope, int i_component);

    double EFT_formfactor(double E, int i_op, int i_isotope, int i_component);

    //EFT response functions
	double F_EFT (double E, int i_isotope, int i_component, int i_op);
	
	//Wrapper functions for each of the different i_fun components
	double F_M(double E, int i_isotope, int i_component);
	double F_Sigma1(double E, int i_isotope, int i_component);
	double F_Sigma2(double E, int i_isotope, int i_component);
	double F_Delta(double E, int i_isotope, int i_component);
	double F_Phi2(double E, int i_isotope, int i_component);
	double F_MPhi2(double E, int i_isotope, int i_component);
	double F_Sigma1Delta(double E, int i_isotope, int i_component);

	double F_SD(double E, int i_isotope, int i_component);

    //Enhancement factors
    double SI_enhancement(int i_isotope);
    double SD_enhancement(int i_isotope);

	//Experimental performance functions
	double efficiency(double E);

    //-------------------------------ADD HERE THE ROUTINE THAT LOADS IN THE DATA???----------

    void print_data(std::string filename);
    void load_data(std::string filename);
    void print_asimov_data(std::string filename);
    void load_asimov_data(std::string filename);

    //Routine for binning the data - ONCE it has been loaded in
    void bin_data();

};


#endif
