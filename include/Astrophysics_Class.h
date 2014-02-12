#ifndef ASTROCLASS_H
#define ASTROCLASS_H

#include <string>

class Astrophysics
{

    public:
        static double rho_x;
        std::string dist_type;           

        //Maxwell-Boltzmann parameters (allowing for multiple component)
        int N_dist;
        double* fraction;
        double* v_lag; //Only allowing scalar v_lag...
        double* v_rms;

        //Lisanti et al. parameters
        double v0;
        int k; //Only using integer spectral index at the moment

        //Parametrization
        int N_vp;
        double* vel_params;  

        double* vel_params_forward;
        double* vel_params_backward;

        //Binned parametrisation
        double* bin_edges;
        double* forward_bin_edges;
        double* backward_bin_edges;

        double initialise_bins(int N, int dir);
        int normalise_bins();
        double rescale_bins(int dir);


        //Poly-exp parametrisation
        double initialise_terms(int N, int dir);
        double normalise_terms(int dir);

        //General parameters
        double v_esc;

        //-------Function prototypes----------
        double load_params();

        double calc_bin_edges(double start, double end, int N_bins, double* edges);

        double (*velocityIntegral) (double,Astrophysics*);

        Astrophysics();
        ~Astrophysics();

   private:
      static const double v_max;

};



#endif
