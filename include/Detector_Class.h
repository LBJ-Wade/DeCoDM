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

    std::vector<double> N_p; //Number of protons
    std::vector<double> N_n; //Number of neutrons
    
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

    //Bin width and number of bins for binned analysis
    double bin_width;
    double tbin_width;
    int N_Ebins;
    int N_tbins;

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




   //Constructor
    Detector(); //Default constructor
    Detector(std::string); //Initialise from an experimental parameter file

    void displayParameters();
    void displayEvents();

    //Nuclear form factors
    double SI_formfactor(double E, int i_isotope);
    double SD_formfactor(double E, int i_isotope, int i_component);

    //Enhancement factors
    double SI_enhancement(int i_isotope);
    double SD_enhancement(int i_isotope);

    //-------------------------------ADD HERE THE ROUTINE THAT LOADS IN THE DATA???----------

    void print_data(std::string filename);
    void load_data(std::string filename);
    void print_asimov_data(std::string filename);
    void load_asimov_data(std::string filename);

    //Routine for binning the data - ONCE it has been loaded in
    void bin_data();

};


#endif
