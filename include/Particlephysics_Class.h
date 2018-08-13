#ifndef THEORYCLASS_H
#define THEORYCLASS_H

#include <iostream>
#include <vector>

class Particlephysics
{
    public:
        double m_x;
        double sigma_SI;
        double sigma_SD;
		
		double lambda_p_D;
		double lambda_n_D;
		
		double lambda_p_Dbar;
		double lambda_n_Dbar;
		
		double sigma_O5;
		double sigma_O7;
		double sigma_O15;
		
		std::vector<double> c_p;
		std::vector<double> c_n;

 	   	int op1;
		int op2;
		int N1;
		int N2;
 
        double delta;

		//Particlephysics();

	//Ratios of proton to neutron couplings for SI and SD interactions
	double r_SI;   //f_p/f_n
	double r_SD;   //a_p/a_n

        //Include ratio f_p/f_n

	//------Function prototypes---------

	void PrintAll();

	Particlephysics()
	{
		m_x = 100;
		sigma_SI = 0;
		sigma_SD = 0;
		lambda_p_D = 0;
		lambda_n_D = 0;
		lambda_p_Dbar = 0;
		lambda_n_Dbar = 0;
		
		sigma_SD = 0;
		sigma_O5 = 0;
		sigma_O7 = 0;
		sigma_O15 = 0;
		c_p.resize(12);
		c_n.resize(12);
		
		op1 = 0;
		op2 = 0;
		N1 = 0;
		N2 = 0;
	}


};

#endif





