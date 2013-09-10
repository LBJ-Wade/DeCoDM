#ifndef EVENT_RATES
#define EVENT_RATES

#include "ParamSet_Class.h"

//-----------Function Prototypes----------

double maxwellRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma);
double diffRate (double v, double theta, double phi, ParamSet params);
double maxwell(double E, void* params);

double binnedRate(double v, double theta, double phi, int N_bins, double* v_egdes, ParamSet params);


double multipoleRadon(double v_q_tmp, int l_tmp, double rate (double,void*), double* params);
double multipoleIntegrand(double v, void* params);

double isotropicRateBinned(double E, void* params);
double forwardRateBinned(double E, void* params);
double backwardRateBinned(double E, void* params);

double isotropicRatePoly(double E, void* params );
double forwardRatePoly(double E, void* params );
double backwardRatePoly(double E, void* params );
double polyf(double v, void* params);


double polyf_DIR(double v, void* params);

double multipoleIntegrand_f0(double v, void* params);
double multipoleIntegrand_f1(double v, void* params);

double forwardIntegrandPoly(double v, void* params);
double backwardIntegrandPoly(double v, void* params);

double forwardIntegrandBinned(double v, void* params);
double backwardIntegrandBinned(double v, void* params);

double BGRate(double E, void* params);

double convolvedRate(double E, void* params);
double preConvolvedRate(double E, void* params);
double spreadRate(double E, void* params);

double Lisanti_f(double v, void* params);
double Lisanti_norm(void* params);
double LisantiRate(double E, void* params);
double LisantiIntegrand(double v, void* params);

double integrator(double v_q, double integrand (double,void*), double* params);

#endif
