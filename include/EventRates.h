#ifndef EVENT_RATES
#define EVENT_RATES

//#include "ParamSet_Class.h"

//-----------Function Prototypes----------

double DMRate(double E, void* params);
double DMRateNR(double E, void* params);
double DMRateDirectional(double E, double theta, double phi, void* params);
double BGRate(double E, void* params);

double convolvedRate(double E, void* params);
double preConvolvedRate(double E, void* params);
double spreadRate(double E, void* params);

double maxwellRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma);
double maxwellModifiedRadon(double v, double theta, double phi, double vlagx, double vlagy, double vlagz, double sigma);
//double diffRate (double v, double theta, double phi, ParamSet params);
//double maxwellRate(double E, void* params);


#endif
