#ifndef PARAMSETCLASS_H
#define PARAMSETCLASS_H

#include "Detector_Class.h"

class ParamSet
{
  public:
    Detector* exptParams;
    double* theoryParams;
  
  //Default constructor
    ParamSet(Detector* expt, double* theory)
	    :exptParams(expt), theoryParams(theory) {};
  
};


#endif
