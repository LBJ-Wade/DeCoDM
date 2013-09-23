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

    double m_det; //Total detector mass in kg
    double exposure; //Effective exposure (t_exp x eff)
    double E_min;
    double E_max;


    //Spin parameters
    double Sn;
    double Sp;
    double zeta;
    double N;
    double alpha;
    double beta;

    double J; //Total Spin

    //Background event levels (in events /keV/kg/day
    double BG_level;

    //Bin data and use binned likelihood analysis
    int USE_BINNED_DATA;

    //Bin width and number of bins for binned analysis
    double bin_width;
    int N_Ebins;

    //Energy resultion (keV)
    double dE;


    std::vector<Event> data; //List of events in detector
    std::vector<double> bin_edges; //List of energy bin edges
    std::vector<int> binned_data; //List of binned data
    std::vector<double> asimov_data; //List of asimov data in each bin
    double No(){return data.size();}




   //Constructor
    Detector(); //Default constructor
    Detector(std::string); //Initialise from an experimental parameter file

    void displayParameters();
    void displayEvents();

    //Nuclear form factors
    double SI_formfactor(double E, int i_isotope);
    double SD_formfactor(double E, int i_isotope);

    //Enhancement factors
    double SI_enhancement(int i_isotope);
    double SD_enhancement();

    //-------------------------------ADD HERE THE ROUTINE THAT LOADS IN THE DATA???----------

    void print_data(std::string filename);
    void load_data(std::string filename);
    void print_asimov_data(std::string filename);
    void load_asimov_data(std::string filename);

    //Routine for binning the data - ONCE it has been loaded in
    void bin_data();

};


#endif
