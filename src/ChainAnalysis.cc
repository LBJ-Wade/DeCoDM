#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "DMUtils.h"

extern "C" { void dsinterface_nevents_( double*, double*, double*, double*, double*, double*, double*);}
extern "C" { void dsinterface_init_( double*, int*, int*);}


//--------Main - Event Generator---------                                                              
int main(int argc, char *argv[])
{
  std::string chainID;
  std::string prefix;
  std::string chainfile;

    //Read in command line arguments                                                                     
    if (argc < 3)
      {
	std::cout << "Not enough arguments: need chain ID (e.g. 'pre067') and chain prefix (e.g. 'DMDD_IC_bin')" << std::endl;
	exit (EXIT_FAILURE);
      }

     chainID = argv[1];
     prefix = argv[2];

     chainfile = "/home/bkavanag/chains/JOB_"+chainID;
     chainfile += "/" + prefix + ".txt";

     //Open the chainfile
     //Read in decay channel parameters...                                                                
     std::cout << "Opening chain file:\t " << chainfile << std::endl;

     std::ifstream chain (chainfile.c_str());
     std::string temp;
     std::vector<double> line;
     double value;
     int N_samps = 0;
     int N_params = 0;     
     int N_length = 0;

     std::vector<double> data;
     std::vector<double> posterior;
     std::vector<double> likelihood;
     std::vector<double> mass;
     std::vector<double> sigsi;
     std::vector<double> sigsd;
     std::vector< std::vector<double> > velparams;

     if (chain.is_open())
       {
	 for (N_samps = 0; std::getline(chain, temp); ++N_samps)
	   {
	     std::stringstream buffer(temp);
             while (buffer >> value) data.push_back(value);
	     //std::istringstream buffer(temp);
	     //std::vector<double> line(std::istream_iterator<double>(buffer),
	     //		      std::istream_iterator<double>());


	     //data.push_back(line);
	   }
	 chain.close();
       }
     else 
       {
           std::cout << "Unable to open chainfile file:\t" << chainfile << std::endl;
           exit(EXIT_FAILURE);
       }
     N_params = (data.size()/N_samps)-2;

     std::cout << "N_samps:\t" << N_samps << std::endl; 
     std::cout << "N_params:\t" << N_params << std::endl;
     //std::cout << "First 2 elements of data:\t" << data[0] << "\t" << data[1] << std::endl;

     int N_vp =N_params-3;

     std::vector<double> vp;
     std::vector<double> vp_max;
     std::vector<double> vp_min;
   
     vp_max.resize(N_vp, 0.0);
     vp_min.resize(N_vp, 1.0);

     //Unpack the chain data
     for (int i = 0; i < N_samps; i++)
       {
	 posterior.push_back(data[0 + (N_params+2)*i]);
         likelihood.push_back(data[1 + (N_params+2)*i]);
          
         mass.push_back(data[2 + (N_params+2)*i]);
         sigsi.push_back(data[3 + (N_params+2)*i]);
         sigsd.push_back(data[4+ (N_params+2)*i]);
         if (N_vp > 0)
	   {
	     vp.clear();
              for (int j = 0; j < N_vp; j++)
	      {
                  vp.push_back(data[5+j + (N_params+2)*i]);
		 
                  if (vp[j] > vp_max[j]) vp_max[j] = vp[j];
                  if (vp[j] < vp_min[j]) vp_min[j] = vp[j];
	      }
              velparams.push_back(vp);
           }
       }
         data.clear();
     
	 for (int i = 0; i < N_samps; i++)
	   {
	     //std::cout << mass[i] << std::endl;
	   }

         int index_bf = std::min_element(likelihood.begin(), likelihood.end()) - likelihood.begin();
         double minlike = likelihood[index_bf];
	 std::cout << "Minimum -2*lnL:\t" << minlike << std::endl; 
	 std::cout << "Bestfit point:" << std::endl;
	 std::cout << "\tm_chi = " << pow(10,mass[index_bf]) << std::endl;
	 std::cout << "\tsigma_SI = " << pow(10,sigsi[index_bf]) <<std::endl;
	 std::cout << "\tsigma_SD = " << pow(10,sigsd[index_bf]) <<std::endl;
         for (int j = 0; j < N_vp; j++)
	   {
	     std::cout << "\tvelparams[" << (j+1) << "] = " << (velparams[index_bf])[j] << std::endl;
	   }


	 std::cout << std::endl;
	 std::cout << "Livepoints limits:" << std::endl;
         for (int j = 0; j < N_vp; j++)
	   {
	     std::cout << "\tvelparams[" << (j+1) << "] in " << vp_min[j] << " --> " << vp_max[j] << std::endl;
	   }

	 //double f[10] = {0.67, 0.311, 0.00237, 0.00188, 0.00878, 0.00193, 0.000733, 0.000798, 0.000550, 0.00142};
           
  //fortfunc_(&ii, &ff);

  return 0;
}
