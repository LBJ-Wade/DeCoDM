#ifndef ASTROCLASS_H
#define ASTROCLASS_H

class Astrophysics
{
	//DO I EVEN NEED ANOTHER CLASS - 'DISTRIBUTION'?

    public:
        double rho_x;
        std::string dist_type;    

	double normalise_bins();
        

        //Maxwell-Boltzmann parameters (allowing for multiple component)
        int N_dist;
        double* fraction;
        double* v_lag; //Only allowing scalar v_lag...
        double* sigma_v;

        //Lisanti et al. parameters
        double v_0;
        int k; //Only using integer spectral index at the moment

        //Parametrization
        int N_vp;
        double* vel_params;  

        //General parameters
        double v_esc;

        //-------Function prototypes----------
        double load_params();
};



#endif
