#include "Neutrinos.h"
#include "Detector_Class.h"
#include <fstream>
#include "ParamSet_Class.h"

int main()
{
    Detector expt = Detector("../test/DDexpt/Experiment1.txt");
	ParamSet parameters(&expt,NULL, NULL);
	LoadFluxTable();
	
	double E1 = 0.1;
	double E2 = 100;
	int N = 1000;
	double dE = (E2 - E1)/(1.0*N);
	
	std::ofstream outputfile;
	outputfile.open("NeutrinoRate.txt");
	
	double Ei = 0;
	for (int i = 0; i < N; i++)
	{
		Ei = E1 + i*dE;
		outputfile << Ei << "\t" << NeutrinoRate(Ei, &parameters) << std::endl;
		
	}
	outputfile.close();
	ClearFluxTable();
	return 1;
	
}