#include "Astrophysics_Class.h"
#include "DMUtils.h"
#include <string.h>
#include <stdlib.h>
//USE A VECTOR FOR V_LAG

#include "gsl/gsl_sf.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

#include "Distributions.h"

const double Astrophysics::v_max = 1000.0;
double Astrophysics::rho_x = 0;


Astrophysics::Astrophysics()
{
    fraction = NULL;
    v_lag = NULL;
    v_rms = NULL;
	v_lag_x = NULL;
	v_lag_y = NULL;
	v_lag_z = NULL;

    vel_params = NULL;

    vel_params_forward = NULL;
    vel_params_backward = NULL;

    bin_edges = NULL;
    forward_bin_edges = NULL;
    backward_bin_edges = NULL;


}

//-------------------------------------------------
double Astrophysics::modulated_rate(double t)
{
    //This method calculates the rate correction factor mod_correction[]
    //associated with the time of each 'super'bin.
    
    
    
}


//-------------------------------------------------
double Astrophysics::load_params()
{
  //Open file for reading in distribution parameters
  std::ifstream file ("dist.txt");
  if (file.is_open())
  {
    //Load in general astrophysical parameters
    rho_x = read_param_double(&file, "rho_x");
    v_esc = read_param_double(&file, "v_esc");

      
    //Load in modulation parameters
    //mod_amplitude = read_param_double(&file, "mod_amplitude");
    //mod_phase = read_param_double(&file, "mod_phase");
    //mod_period = read_param_double(&file, "mod_period");
      
    dist_type = read_param_string(&file, "dist_type");
    //std::cout << "Using " << dist_type << " type distribution..." << std::endl;

    //Read in distribution-specific parameters
    if (dist_type == "maxwell")
    {
      char numstr[21]; // enough to hold all numbers up to 64-bits

      //Read in number of distributions and initialise arrays
      N_dist = read_param_int(&file, "N_dist");
      fraction = new double[N_dist];
      v_lag = new double[N_dist];
      v_rms = new double[N_dist];

      v_lag_x = new double[N_dist];
	  v_lag_y = new double[N_dist];
	  v_lag_z = new double[N_dist];

      double v_lag_vec[3];

      for (int i = 0; i < N_dist; i++)
      {
		  sprintf(numstr, "%d", i+1);
		  fraction[i] = read_param_double(&file, "fraction"+std::string(numstr));
		  read_param_vector(&file, "v_lag"+std::string(numstr),v_lag_vec, 3);
		  v_rms[i] = read_param_double(&file, "v_rms" + std::string(numstr));

        //Calculate length of v_lag
        v_lag[i] = sqrt(v_lag_vec[0]*v_lag_vec[0] + v_lag_vec[1]*v_lag_vec[1] + v_lag_vec[2]*v_lag_vec[2]);
		v_lag_x[i] = v_lag_vec[0];
		v_lag_y[i] = v_lag_vec[1];
		v_lag_z[i] = v_lag_vec[2];
      }
      velocityIntegral = &velInt_maxwell;
	  modifiedVelocityIntegral = &velInt_maxwell_modified;
	  //std::cout << "In Astrophysics_Class.cc: " << std::endl;
	  //std::cout << "\t\t Using velInt_maxwell_modified..." << std::endl;
    }
    else if (dist_type == "lisanti")
    {
      //Read in parameter values
      v0 = read_param_double(&file, "v0");
      k = read_param_int(&file, "k");

      velocityIntegral = &velInt_Lisanti;
    }
    else if (dist_type == "param")
    {
      int mode = read_param_int(&file, "mode");
      if (mode == 1)
      {
	int N = read_param_int(&file, "N_vp");
	initialise_bins(N, 0);
	read_param_vector(&file, "vel_params", vel_params, N_vp);
	rescale_bins(0);

	velocityIntegral = &velInt_isotropicBinned;
      }
      else
      {
	std::cout << "Parametrisation mode '" <<  mode << "' is not valid. Exiting..." << std::endl;
	exit (EXIT_FAILURE);
      }

    }
    else
    {
     std::cout << "dist_type '" <<  dist_type << "' is not valid. Exiting..." << std::endl;
     exit (EXIT_FAILURE);
    }
  }
  else
  {
    std::cout << "Unable to open distribution parameter file:\t'" << "dist.txt" << "'" << std::endl;
    exit (EXIT_FAILURE);
  }
  file.close();
}




//-----------------------------------------------------------------
//-------------Binned parametrisation------------------------------
//-----------------------------------------------------------------

double Astrophysics::initialise_bins(int N, int dir)
{
  N_vp = N;

 if (dir == 0)
 {
   //Initialise arrays
    vel_params = new double[N_vp];
    bin_edges = new double[N_vp +1];
    //Calculate bin edges
    calc_bin_edges(0, v_max, N_vp, bin_edges);
 }
 else
 {
   vel_params_forward = new double[N_vp];
   vel_params_backward = new double[N_vp];
   forward_bin_edges = new double[N_vp +1];
   backward_bin_edges = new double[N_vp +1];
   //Calculate bin edges
   calc_bin_edges(0, v_max, N_vp, forward_bin_edges);
   calc_bin_edges(0, v_max, N_vp, backward_bin_edges);

 }
}

//-----------------------------------------WHAT ABOUT DIRECTIONAL BINS...?
int Astrophysics::normalise_bins()
{
   double norm = 0;

   for (int i = 1; i < N_vp; i++)
   {
       norm += vel_params[i];
   }
   if (norm > 1)
   {
       return -1;
   }
   else
   {
       vel_params[0] = 1.0-norm;
       return 1;
   }
}

double Astrophysics::calc_bin_edges(double start, double end, int N_bins, double* edges)
{

    double bin_width = (end - start)/N_bins;

    for (int i = 0; i < N_bins+1; i++)
    {
      edges[i] = start + i*bin_width;
    }
}

double Astrophysics::rescale_bins(int dir)
{
   if (dir == 0)
   {
       //Calculate normalised bin values
       for (int i = 0; i < N_vp; i++)
       {
           vel_params[i] = 3*vel_params[i]/(pow(bin_edges[i+1],3) - pow(bin_edges[i],3));
       }
   }
   else
   {
       //Calculate normalised bin values
       for (int i = 0; i < N_vp; i++)
       {
           vel_params_forward[i] = 3*vel_params_forward[i]/(pow(forward_bin_edges[i+1],3) - pow(forward_bin_edges[i],3));
           vel_params_backward[i] = 3*vel_params_backward[i]/(pow(backward_bin_edges[i+1],3) - pow(backward_bin_edges[i],3));
       }
   }
}




//-----------------------------------------------------------------
//-------------Polynomial-Exponential parametrisation--------------
//-----------------------------------------------------------------

double Astrophysics::initialise_terms(int N, int dir)
{
   N_vp = N;
   if (dir == 0)
   {
      vel_params = new double[N_vp];
   }
   else
   {
	   std::vector<double> vp(N_vp);
	   for (int i = 0; i < N_ang; i++)
	   {
		   vel_params_ang.push_back(vp);
		
	   }
   }
}

double Astrophysics::normalise_terms(int dir, double* norms)
{
  //Declare gsl workspace (1000 subintervals)
  gsl_integration_workspace * workspace
	      = gsl_integration_workspace_alloc (3000);

  //Declare gsl function to be integrated
  gsl_function F;

  double norm = 0;
  double norm2 = 0;
  if (dir == 0)
  {
      F.function = &polyf;
	  F.params = NULL;

	  double error;
	  int status = gsl_integration_qag(&F,0,v_max, 0, 1e-6, 3000,6,workspace, &norm, &error);

	  if (status ==  GSL_EROUND)
	  {
	    std::cout << "GSL rounding error!" << std::endl;
	    std::cout << norm << std::endl;
	  }

      vel_params[0] += log(norm);

      //Free workspace
      gsl_integration_workspace_free (workspace);
  }
  else
  {
      //------------------------------------------------
      //std::cout  << "Normalisation not working for non-directional distributions - need to distinguish between forward and back...Also need to actually adjust - possibly use a ratio parameter" << std::endl;
       //exit (EXIT_FAILURE);
      //-------------------------------------------------

      double error;
      int status;
	  
	  double totalnorm = 0;	
	  double totalnorm2 = 0;
      F.function = &polyf_angnorm;
      F.params = this;

      
      double fudge = 1.0;
	  for (int k = 0; k < N_ang; k++)
	  {

	    double v_cut = 0;
		 norm = 0;
		j_bin = k;
	  	status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,workspace, &norm, &error);
		status = gsl_integration_qag(&F,v_cut,1000, 0, 1e-6, 3000, 6, workspace, &norm2, &error);
        if (status ==  GSL_EROUND)
        {
  	    std::cout << "GSL rounding error!" << std::endl;
  	    std::cout << norm << std::endl;
        }
	
	fudge = 1.0;

	//Fudged for folded...
	if (k == 0) fudge = 1.0;
	if (k == 2) fudge = 0.0;
	norms[k] = fudge*norm2;
                totalnorm2 += fudge*norm2;
		totalnorm += fudge*norm;
	  }
	  
	  //This is a strangely fixed normalisation!
	  for (int k = 0; k < N_ang; k++)
	  {
		  //vel_params_ang[k][0] += log(N_ang*norm*(cos(PI*(k)/N_ang) - cos(PI*(k+1.0)/N_ang)));
		  norms[k] = norms[k]/totalnorm2;
		  vel_params_ang[k][0] += log(totalnorm);
	  }
	  
      //Free workspace
      gsl_integration_workspace_free (workspace);
  }
}


/*
double velocityIntegral(double v, void* params)
{
    if (dist_type == "maxwell")
    {
        return velInt_maxwell(v, params);
    }
    else if (dist_type == "lisanti")
    {
        return velInt_Lisanti(v, params);
    }
    else if (dist_type == "binned")
    {
        return velInt_isotropicBinned(v, params);
    }
    else if (dist_type == "polyexp")
    {
        return velInt_isotropicPoly(v, params);
    }
}*/


// Destructor
Astrophysics::~Astrophysics()
{
    //Free up memory
    delete[] fraction;
    delete[] v_lag;
    delete[] v_rms;

    delete[] bin_edges;
    delete[] forward_bin_edges;
    delete[] backward_bin_edges;
    delete[] vel_params;
    delete[] vel_params_forward;
    delete[] vel_params_backward;

}
