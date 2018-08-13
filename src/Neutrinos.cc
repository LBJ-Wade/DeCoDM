#include <math.h>
#include <iostream>

#include "Detector_Class.h"
#include "ParamSet_Class.h"
//#include "EventRates.h"
#include "DMUtils.h"
#include "Neutrinos.h"

//#include "Astrophysics_Class.h"
//#include "Particlephysics_Class.h"

//#include "gsl/gsl_sf.h"
//#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


double E_R_temp;
int NeutrinoID;
int i_isotope;

//Tabulated neutrino fluxes
//double* E_range;
//double* N_flux;

//End points of tabulated energies
//double E_start;
//double E_end;
//int N_lines;

//Interpolation machinery
gsl_interp_accel *acc;
gsl_spline *spline;


//---------------------------------------------------------------------------
//------------Neutrino background rates--------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------
//Background rate in events /keV/kg/day
double NeutrinoRate(double E, void* params)
{
	//return 1e-15;
	
	E_R_temp = E;
	
    //Declare gsl workspace (1000 subintervals)
    gsl_integration_workspace * workspace= gsl_integration_workspace_alloc (5000);

    //switch off default error handler, store old error handler in
    //old_handler
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

    //Declare gsl function to be integrated
    gsl_function F;


    F.params = params;

    Detector* expt = ((ParamSet*)params)->exptParams;

    double result, error;
    int status;
     
	F.function = NeutrinoIntegrand;

    double IntegralResult = 0;

	for (int j = 0; j < expt->N_isotopes; j++)
	{
		i_isotope = j;
		//Need to check for E_min_v below lower limit!
        double Ev1 = E_min_v(E, expt->m_n[j]);
		//std::cout << Ev1 << std::endl;
		//std::cout << Ev1 << std::endl;
		//if (Ev1 < 100) Ev1 = 101;
	    double Ev2 = 1e6; //100 MeV = 1e5 keV

		//Calculate weak nuclear hypercharge
		double Q = expt->N_n[j] - (1.0-4.0*0.2397)*expt->N_p[j];	

        status = gsl_integration_qag(&F,Ev1,Ev2, 0, 1e-6, 5000,6, workspace, &result, &error);
	    IntegralResult += expt->frac_n[j]*result*Q*Q*expt->SI_formfactor(E,i_isotope);
		//std::cout << Q << "\t" << Ev1 << "\t" << result << std::endl;
    }			 
				
	//Fermi coupling constant (in keV^-2)
	double G_F = 1.16637e-17;
	
	return 1.88e22*(G_F*G_F)*IntegralResult/(4.0*PI);
							   
}

//-----------------------------------------------------
//Calculate integrand that appears in the integral over neutrino energies
double NeutrinoIntegrand(double E_v, void* params)
{
	Detector* expt = ((ParamSet*)params)->exptParams;
	
	
	//Calculate nuclear mass in keV
	double m_n = expt->m_n[i_isotope]*(931.5e-3)*1e6;
	//Need to get the rescaling of the units correct!
	
	//Need to put in here something about the upper recoil energy!!!
	
	return (1.0-m_n*E_R_temp/(2.0*E_v*E_v))*NeutrinoFlux(E_v);
}

//------------------------------------------------------
//Minimum neutrino energy required for nuclear recoil (all in keV)
double E_min_v(double E, double m_n)
{
	//Nuclear mass in keV
	m_n *= (931.5e-3)*1e6;
	
	return sqrt(0.5*E*m_n);
}

//-----------------------------------------------------
//Differential neutrino flux
double NeutrinoFlux(double E_v)
{
	//May need a factor of 1000 or something here to go from 1/MeV to 1/keV!
	
	//May need to check - am I integrating in log-space or lin-space?
	
    //Need to do checks for 'out of range' energy values?
	
	/*
	if ((E_v > 0.99e4) && (E_v < 1.01e4))
	{
		return 1;
	}
	else
	{
		return 0;
	}
	*/
	return 1e-3*gsl_spline_eval (spline, E_v, acc); 
	
}



//-----------------------------------------------------
//Load neutrino flux table from file
double LoadFluxTable()
{
	//std::cout << "Need a flag to specify whether the flux table has been loaded..." << std::endl;
	std::string line;
    std::ifstream file ("data/NeutrinoSpectrum.txt");
    if (file.is_open())
    {
	    //Calculate number of lines
	    int N_lines = 0;
		double temp = 0;
	    while (std::getline(file, line))
	   	{
	           ++N_lines;
	   	}
	
	    //Reset to start of file
	    file.clear();
	    file.seekg(0, std::ios::beg);
	
	    //NB:remember to pull...? Check where the stuff is...
        double* E_range;
		double* N_flux;
			
		E_range = (double*)malloc(N_lines*sizeof(double));
		N_flux = (double*)malloc(N_lines*sizeof(double));
		
		std::string val;
	    for (int i = 0; i < N_lines; i++)
	    {
			file >> val;
			E_range[i] = 1e3*atof(val.c_str()); //In keV
			file >> val;
			N_flux[i] = atof(val.c_str());
			//std:: cout << E_range[i] << '\t' << N_flux[i] << std::endl;
			//Do I need to check for new lines???
	    }
		
		acc = gsl_interp_accel_alloc ();
		spline = gsl_spline_alloc (gsl_interp_cspline, N_lines);

	    gsl_spline_init (spline, E_range, N_flux, N_lines);
		
		free(E_range);
		free(N_flux);
	    file.close();
		return 1;
		
    }
	else
	{
		std::cout << "Neutrino file - 'data/NeutrinoSpectrum.txt' - not found..." << std::endl;
		return 0;
	}
}

//----------------------------------------------------
//Clear memory
double ClearFluxTable()
{
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
	
}

