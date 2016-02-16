#ifndef GENERATEEVENTS_H
#define GENERATEEVENTS_H

#include "Detector_Class.h"
#include <vector>
#include <string.h>

//--------Function Prototypes-----------

void generateEvents(Detector* expt,double m_x, double sigma_SI, double sigma_SD);
void generateEvents_Ne(Detector* expt, double m_x, int No_target, int op);
double calcNe(Detector* expt, double m_x, double sigma_SI, double sigma_SD);
double calcNe_NR(Detector* expt, double m_x, int op1, int op2, int N1, int N2);
void printEvents(std::vector<Event> data, std::string filename);

double generateBGEvents(Detector* expt);
void generateBGEvents_Asimov(Detector* expt);

double generateNeutrinoEvents(Detector* expt);
void generateNeutrinoEvents_Asimov(Detector* expt);

void calcRotationMatrix(double* rot_matrix, double* v_lag);
void rotateEvent(double* theta, double* phi, double* rot_matrix);

void loadExperiments();
void initialiseRNG();
void clearRNG();
void generateAllEvents(double m_x, double sigma_SI, double sigma_SD);
void generateAllEvents_Ne(double m_x, int No_target, int op);
void printAllData();
	 
#endif
