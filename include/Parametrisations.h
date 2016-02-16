#ifndef PARAMETRISATIONS_H
#define PARAMETRISATIONS_H


double f0_binned(double v,double* g,double* v_edges, int N_bins);
double f1_binned(double v,double* g,double* v_edges, int N_bins);
double f2_binned(double v,double* g,double* v_edges, int N_bins);
double eta_binned(double v,double* g,double* v_edges, int N_bins);

//Directional discretisation
double J(double beta, double x, double y);
double fintegrand(double v, void * params);
double velInt_DRT(double vmin, Astrophysics* astro);
#endif
