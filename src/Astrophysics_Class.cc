#include "Astrophysics_Class.h"
#include "DMUtils.h"

//USE A VECTOR FOR V_LAG

double Astrophysics::load_params()
{
  //Open file for reading in distribution parameters
  std::ifstream file ("dist.txt");
  if (file.is_open())
  {
    rho_x = read_param_double(&file, "rho_x");
    v_esc = read_param_double(&file, "v_esc");

    dist_type = read_param_string(&file, "dist_type");
    std::cout << "Using " << dist_type << " type distribution..." << std::endl;

    if (dist_type == "maxwell")
    {
      //Read in number of distributions and initialise arrays
      N_dist = read_param_int(&file, "N_dist");
      fraction = new double[N_dist];      
      v_lag = new double[N_dist];
      sigma_v = new double[N_dist];

      for (int i = 0; i < N_dist; i++)
      {
	sprintf(numstr, "%d", i+1);
	fraction[i] = read_param_double(&file, "fraction"+std::string(numstr));
	read_param_vector(&file, "v_lag"+std::string(numstr),v_lag);
	sigma_v[i] = read_param_double(&file, "sigma_v" + std::string(numstr));

	//std::cout << sigma_v << '\t' << v_lag[2] << std::endl;
	Ne += generateMaxwellEvents(expt, m_x, fraction*sigma_SI, fraction*sigma_SD, v_lag, sigma_v, v_esc);
	generateMaxwellEvents_Asimov(expt, m_x, fraction*sigma_SI, fraction*sigma_SD, v_lag, sigma_v, v_esc);

      }
    }
}

// Destructor
Astrophysics::~Astrophysics() 
{
    //Free up memory
    delete[] fraction; 
    delete[] v_lag;
    delete[] sigma_v;

}
