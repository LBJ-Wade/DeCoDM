#ifndef DISTCLASS_H
#define DISTCLASS_H

class Distribution
{
	//DO I EVEN NEED ANOTHER CLASS - 'DISTRIBUTION'?

    public:
        double rho_x;
        
       
        int N_vp;
        double* vel_params;        

        //Distribution parameters
        double* frac;

        //Maxwell-Boltzmann parameters (allowing for multiple component)

        //Lisanti et al. parameters
        double v_0;
        int k; //Only using integer spectral index at the moment

        //General parameters
        double v_esc;

        //-------Function prototypes----------
        double dist_function();
        double EventRate();
};



#endif
