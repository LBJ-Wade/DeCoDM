#ifndef DM_UTILS_H
#define DM_UTILS_H

#include <iostream>
#include <string>
#include <math.h>

#include "ParamSet_Class.h"

#ifndef PI
  #define PI 3.14159265358
#endif

//------------Function Prototypes---------------

double reduced_m(double m_n, double m_x);
double reduced_m_GeV(double m_n, double m_x);
double rate_prefactor(double m_n, double m_x, double sigma, double rho);

double v_min(double E, double m_n, double m_x);
double p_min(double E, double m_n);

double u_max(double v, double m_n, double m_x);

//double* calcBins(double start, double end, int N_bins);

double v_min_inverse(double v, double m_n, double m_x);

double N_expected(double rate (double,void*), ParamSet parameters);
double N_expected(double rate (double,void*), ParamSet parameters, double E1, double E2);

double sphericalHarmonic(int l, int m, double theta, double phi);

double logfactNo(double No);

double PoissonLike(Detector* expt, ParamSet parameters, double signal_rate (double,void*), double No, double E1, double E2);
double signal_fraction(Detector* expt, ParamSet parameters, double signal_rate (double, void*));

//Routines for reading in parameters from files
int read_param_int(std::ifstream* file, std::string param_name);
double read_param_double(std::ifstream* file, std::string param_name);
double read_param_vector(std::ifstream* file, std::string param_name, double* output, int length);
std::string read_param_string(std::ifstream* file, std::string param_name);

void setCurrentRate( double rate(double, void*));
void setCurrentVelInt( double VelInt(double, void*));

int load_params(std::string filename);

double ChebyshevP(int order, double x);

int nCr(int n, int r);

double ChebyshevP(int order, double x);

//------------Global Variables

extern std::string expt_folder;
extern std::string events_folder;
extern int N_expt; //Number of experiments
extern int vmode;   //Speed parametrization mode (see guide...)
extern int USE_SD; //Use spin-dependent interactions
extern int USE_SI; //Use spin-dependent interactions
extern int DIR;	   //Use directional information on events
//extern int N_terms; //How many terms are used in the parametrisation
extern int USE_FLOAT_BG; //Use a floating background level in each experiment
extern int USE_VARY_FF; //Use varying SD Form Factor in each experiment
extern int USE_ASIMOV_DATA; //Analyse using asimov data

//Define globally stored limits for Ne integration
extern double E_a;
extern double E_b;

//Variables needed to make the integration simpler (so they don't need storing as pointers)
extern double currentE;
extern double (*currentRate) (double, void*);


//-----------Inline functions...



#endif
