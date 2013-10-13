#ifndef PARAMSETCLASS_H
#define PARAMSETCLASS_H

#include "Detector_Class.h"
#include "Particlephysics_Class.h"
#include "Astrophysics_Class.h"

class ParamSet
{
  public:
    Detector* exptParams;
    Particlephysics* theoryParams;
    Astrophysics* astroParams;
  
  //Default constructor
    ParamSet(Detector* expt, Particlephysics* theory, Astrophysics* astro)
	    :exptParams(expt), theoryParams(theory), astroParams(astro) {};
  
};


#endif
