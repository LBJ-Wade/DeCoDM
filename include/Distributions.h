#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "Astrophysics_Class.h"

    double multipoleRadon(double v_q_tmp, int l_tmp, double rate (double,void*), double* params, int N);

	double velInt_maxwell(double v, Astrophysics* astro);
    double velInt_maxwell_modified(double v, Astrophysics* astro);
  
    double velInt_Lisanti(double v, Astrophysics* astro);

	double velInt_isotropicBinned(double v, Astrophysics* astro);
	double velInt_forwardBinned(double v, Astrophysics* astro);
	double velInt_backwardBinned(double v, Astrophysics* astro);

	double velInt_isotropicPoly(double v,  Astrophysics* astro );
	double velInt_forwardPoly(double v,  Astrophysics* astro );
	double velInt_backwardPoly(double v , Astrophysics* astro );
	

	//double polyf_DIR(double v, void* params);

	//double multipoleIntegrand_f0(double v, void* params);
	//double multipoleIntegrand_f1(double v, void* params);

	double forwardIntegrandPoly(double v, void* params);
	double backwardIntegrandPoly(double v, void* params);

	double forwardIntegrandBinned(double v, void* params);
	double backwardIntegrandBinned(double v, void* params);

        double LisantiIntegrand(double v, void* params);
        double multipoleIntegrand(double v, void* params);

	double Lisanti_f(double v, void* params);
        double Lisanti_norm(void* params);
		
		double velInt_polytotal(double v, Astrophysics* astro);
		
        double polyf(double v, void* params);
		double polyf_ang(double v, Astrophysics* astro, int k);
		double polyf_total(double v, void* params);
		double polyfintegrand_total(double v, void* params);

	double integrator(double v_q, double integrand (double,void*), Astrophysics* astro);

#endif
