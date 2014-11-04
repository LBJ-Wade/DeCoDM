#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "shared.h"
#include "gsl/gsl_integration.h"


void dsinterface_nevents2_( double*, double*, double*, double*, int*, int*, double* ); 
void dsinterface_init_( double*, int*, int*);

double polyf(double v, void* params);

int Feedback;
double *ScalesMin;
double *ScalesMax;

double read_param_double(char *filename, char *param_name);
int read_param_int(char *filename, char *param_name);


//BJK
int vmode;
int FLOAT_BG;
int N_expt;
int USE_VARY_FF;

int USE_ICECUBE;
int USE_DD;

int USE_SI;
int USE_SD;

int N_terms;

double polyf(double v, void* params);
double polyf1(double v, void* params);

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double* INSlogZ, double *logZerr, void *context)
{
  // convert the 2D Fortran arrays to C arrays

  // the posterior distribution
  // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

  int i, j;
  double postdist[*nSamples][*nPar + 2];
  for( i = 0; i < *nPar + 2; i++ )
    for( j = 0; j < *nSamples; j++ )
      postdist[j][i] = posterior[0][i * (*nSamples) + j];

  // last set of live points
  // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

  double pLivePts[*nlive][*nPar + 1];
  for( i = 0; i < *nPar + 1; i++ )
    for( j = 0; j < *nlive; j++ )
      pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}



void LogLike(double *Cube, int *ndim, int *npars, double *Ne)
{
  int i = 0;

  static int count = 0;

  double norm = 0;


  //May need to add to the offset here if we're using the varying form factors...
  int OFFSET = USE_SI+USE_SD+N_expt*FLOAT_BG; 

  if (N_expt == 1) OFFSET += 6*USE_VARY_FF;
  if (N_expt == 2) OFFSET += 6*USE_VARY_FF;
  if (N_expt == 3) OFFSET += 9*USE_VARY_FF;  

  double *Ne_DD;
  Ne_DD = (double*)malloc(N_expt*sizeof(double));
  double Ne_IC = 0;
  
  if (vmode==1)
    {
      
      for (i = (1+OFFSET); i < *npars; i++)
	{
	  norm += Cube[i];
	}
      
      //printf("%f\n", norm);
      if (norm > 1)
	{
	  //*lnew = -1e30;
          //printf("Rejected...\n");
          return;
	}
      
    }

  double GetLogLike = 0.0;
  double GetLogLike_DD = 0.0;

  
  for (i = 0; i < *npars; i++)
    {
      //printf("%f\n", Cube[i]);
    }
  //printf("\n");
  

  for (i = 0; i < *ndim; i++)
    {
      //Cube[i] = ScalesMin[i] + (ScalesMax[i] - ScalesMin[i]) * Cube[i];
      //printf("%f\t", Cube[i]);
    }
  // printf("\n");
  double mchi = pow(10.0, Cube[0]);
  double sigSI = 0;
  if (USE_SI) sigSI = pow(10.0, Cube[1]);
  double sigSD = 0;
  //if ((*ndim > 2)&&(USE_SD)) sigSD = pow(10.0, Cube[2]);
  if ((USE_SI)&&(USE_SD))
    {
      sigSD = pow(10.0, Cube[2]);
    }
  else if (USE_SD)
    {
      sigSD = pow(10.0, Cube[1]);
    }

  //Start BJK
  //Pack up velocity parameters - VelParams[0] will be calculated here and is set by normalisation

  double* VelParams;
  double* DDParams;

  double* VelParams2;

  int N_vp = 0;
  int ndim_new = 0;
  
  gsl_integration_workspace* workspace;
  
  gsl_function F;
  double error;
  int status;

  int N_approx_bins = 1000;


     switch (vmode)
     {
        case 0:
          //Cube[0] = log10(30);
          //Cube[1] = log10(1e-45);
          //Cube[2] = log10(2e-40);
	  if (USE_DD) geteventnumbers_(Cube, ndim, Ne_DD);
           break;   

        case 1:	//Normalisation for binned method
      
	
           N_vp = *npars-OFFSET; //Set number of velocity parameters
           VelParams = (double*)malloc((N_vp)* sizeof(double));

           VelParams[0] = 0;
           for (i = 1; i < (N_vp); i++)
	   {
	      VelParams[i] = Cube[i+OFFSET];
	   }

           //for (i = 0; i < N_vp; i++)
	     //{
	     //VelParams[i] = -log(VelParams[i]);
             //norm += VelParams[i];
	     //}
           
           VelParams[0] = 1 - norm;
           for (i = 0; i < N_vp; i++)
           {
             //printf("%f\t", VelParams[i]);
	     //VelParams[i] /= norm;
           }
	  
           //printf("\n\n");
       
           //Check that VelParams[0] > 0
           //if (VelParams[0] < 0)
	     //{
	     //free(VelParams);
	       //*lnew = -1e30;
		// return;
	       //}



           if (USE_DD) geteventnumbers_(Cube, npars, Ne_DD);
           break;
              
        case 3:	//Normalisation for exponential-polynomial method
 
           ndim_new = *ndim + 1;

           N_vp = *ndim-OFFSET;
           VelParams = (double*)malloc((N_vp)*sizeof(double));
      
           VelParams[0] = 0;
           for (i = 1; i < (N_vp); i++)
           {
              VelParams[i] = Cube[i+OFFSET];
           }
     

           //Declare gsl workspace (1000 subintervals)
           workspace = gsl_integration_workspace_alloc (3000);

           N_terms = N_vp;


           //Declare gsl function to be integrated
           F.function = &polyf;
           F.params = VelParams;

           status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,
				workspace, &norm, &error);

           VelParams[0] = log(norm);


      
           DDParams = (double*)malloc((*ndim + 1)*sizeof(double));
           for (i = 0; i < (1+OFFSET); i++)
	   {
	      DDParams[i] = Cube[i];
	   }

	   //Free workspace
           gsl_integration_workspace_free (workspace);

           /*
           VelParams[0] = 20.7083;
           VelParams[1] = 7.5033;
           VelParams[2] = 1.1354;
           VelParams[3] = 0.6176;
           VelParams[4] = -0.2069;
           VelParams[5] = -0.0490;
           VelParams[6] = 0.1571;
           VelParams[7] = -0.1496;
           VelParams[8] = 0.0808;
           VelParams[9] = -0.0057;
           */

           for (i = 0; i < N_vp; i++)
	     {
	       DDParams[i+1+OFFSET] = VelParams[i];
	     }

           /*
           mchi = 50;
           sigSI = 1e-43;
           sigSD = 2e-40;
           DDParams[0] = log10(30);
           DDParams[1] = log10(1e-45);
           DDParams[2] = log10(2e-40);
           VelParams[0] = 20.7083;
           VelParams[1] = 7.5033;
           VelParams[2] = 1.1354;
           VelParams[3] = 0.6176;
           VelParams[4] = -0.2069;
           VelParams[5] = -0.0490;
           VelParams[6] = 0.1571;
           VelParams[7] = -0.1496;
           VelParams[8] = 0.0808;
           VelParams[9] = -0.0057;
	   */

           
           
           if (USE_DD) geteventnumbers_(DDParams, &ndim_new, Ne_DD);
           break;

     case 4: //Normalisation for approximate exponential-polynomial method

       ndim_new = 1 + OFFSET + N_approx_bins;
                   

       N_vp = *ndim-OFFSET;
       VelParams = (double*)malloc((N_vp)*sizeof(double));

       VelParams[0] = 0;
       for (i = 1; i < (N_vp); i++)
	 {
	   VelParams[i] = Cube[i+OFFSET];
	 }


       //Declare gsl workspace (1000 subintervals)
       workspace = gsl_integration_workspace_alloc (3000);

       N_terms = N_vp;

       //Declare gsl function to be integrated
       F.function = &polyf; //-----------Error here - need to include other lot of code...
       F.params = VelParams;

       status = gsl_integration_qag(&F,0,1000, 0, 1e-6, 3000,6,
					workspace, &norm, &error);

       VelParams[0] = log(norm);

       //printf("%f\n", VelParams[0]); 

       //DDParams = (double*)malloc((*ndim + 1)*sizeof(double));
       DDParams = (double*)malloc((1+OFFSET+N_approx_bins)*sizeof(double));
       for (i = 0; i < (1+OFFSET); i++)
	 {
	   DDParams[i] = Cube[i];
	 }
       
       /*
       mchi = 30;
       sigSI = 1e-45;
       sigSD = 2e-40;
       DDParams[0] = log10(30);
       DDParams[1] = log10(1e-45);
       DDParams[2] = log10(2e-40);
       VelParams[0] = 20.7083;
       VelParams[1] = 7.5033;
       VelParams[2] = 1.1354;
       VelParams[3] = 0.6176;
       VelParams[4] = -0.2069;
       VelParams[5] = -0.0490;
       VelParams[6] = 0.1571;
       VelParams[7] = -0.1496;
       VelParams[8] = 0.0808;
       VelParams[9] = -0.0057;
       */

       double tmp = 0;

       double dv = 1000.0/N_approx_bins;
       
       F.function = &polyf;

       // double tmpsum = 0;

       VelParams2 = (double*)malloc(N_approx_bins*sizeof(double));
       for (i = 0; i < N_approx_bins; i++)
	 {
	   status = gsl_integration_qag(&F, i*dv, (i+1)*dv, 0 , 1e-6, 3000, 6,
					workspace, &tmp, &error);
	   VelParams2[i] = 3.0*tmp/(pow((i+1)*dv, 3) - pow(i*dv, 3));
           //VelParams2[i] = tmp/dv;
           //tmpsum += VelParams2[i]*(pow((i+1)*dv, 3) - pow(i*dv, 3))/3.0;
           //tmpsum += VelParams2[i];	   
           //VelParams2[i] = polyf((i+0.5)*10.0, VelParams);
	 }

       //printf("(Sum of g_i) - 1 = %f\n", tmpsum-1.0);

       //Free workspace
       gsl_integration_workspace_free (workspace);


       for (i = 0; i < N_approx_bins; i++)
	 {
	   DDParams[i+1+OFFSET] = VelParams2[i];
           //DDParams[i+1+OFFSET] = VelParams[i];
           //printf("%d\t %f\n",i, DDParams[i+1+OFFSET]);
	 }


       if (USE_DD) geteventnumbers_(DDParams, &ndim_new, Ne_DD);
       break;

      }
  



  if (USE_ICECUBE)
    {
      if (vmode == 4)
	{
          N_vp = N_approx_bins;
	  //int new_vmode = 3;
          //printf("%d\n", N_vp);
	  dsinterface_nevents2_(&mchi, &sigSI, &sigSD, VelParams2, &N_vp, &vmode, &Ne_IC);
	}
      else
	{
	  dsinterface_nevents2_(&mchi, &sigSI, &sigSD, VelParams, &N_vp, &vmode, &Ne_IC);
	}
    }
 
  //Calculate likelihood contribution from Direct Detection experiments...

  //loglikelihood_(Cube, ndim, &GetLogLike_DD);

  //printf("NB: Using reweighting of bins...\n");

  
   if(vmode > 1) free(VelParams);
   if(vmode == 3) 
     {
        free(DDParams);    
     }
   if (vmode == 4)
     {
       free(DDParams);
       free(VelParams2);
     }
   //printf("Direct detection ln(L): %f\n", -GetLogLike_DD);
   //printf("IceCube ln(L): %f\n", GetLogLike);
   //printf(" %lf %lf %lf\n\n", mchi, log10(sigSI), log10(sigSD));


   //printf("ln(L): %f\n", GetLogLike - GetLogLike_DD);

   //count++;
   //printf("Sampled: %d\n", count);

   //if (count == 1000)
   //  {
   //    exit(EXIT_FAILURE);
   //  }

   double N_tot = 0.0;
   if (N_expt == 1) N_tot = 20.0;
   if (N_expt == 2) N_tot = 55.0;
   if (N_expt == 3) N_tot = 100.0;
   if (USE_ICECUBE == 1) N_tot += 1.0;


   int Ni;
   for (Ni = 0; Ni < N_expt; Ni++)
   {
       Ne[Ni] = Ne_DD[Ni];
       //Ne[Ni] = 0;
   } 
   Ne[N_expt] = Ne_IC;

   free(Ne_DD);


   // *lnew = (N_tot/(N_expt+USE_ICECUBE))*(GetLogLike - GetLogLike_DD);
}



int main(int argc, char *argv[])
{
  int i;
  double logZero = -1e20;
  void *context = 0;

  char root[100];
  sprintf(root, "%s", argv[2]);
  char input_file_name[100];
  sprintf(input_file_name, "%s%s", argv[2], ".ini");
  printf("ini file is %s\n", input_file_name);

  char chain_ID[100];
  sprintf(chain_ID, "%s", argv[1]);
  char chain_file_name[100];
  sprintf(chain_file_name, "%s%s%s%s%s" , "/home/bkavanag/chains/JOB_", argv[1], "/", argv[2], ".txt");
  printf("Chainfile is %s\n", chain_file_name);  

  char chain_file_name_output[100];
  sprintf(chain_file_name_output, "%s%s%s%s%s" , "/home/bkavanag/chains/JOB_", argv[1], "/", argv[2], "2.txt");

  FILE *input_file = fopen(input_file_name,"r");

  fscanf(input_file, "%d", &Feedback);
  printf("Feedback is %d\n", Feedback);

  int fb;
  fscanf(input_file, "%d", &fb);
  printf("fb is %d\n", fb);

  int IS;
  fscanf(input_file, "%d", &IS);
  printf("IS is %d\n", IS);

  int mmodal;
  fscanf(input_file, "%d", &mmodal);
  printf("mmodal is %d\n", mmodal);

  int ceff;
  fscanf(input_file, "%d", &ceff);
  printf("ceff is %d\n", ceff);

  int nlive;
  fscanf(input_file, "%d", &nlive);
  printf("nlive is %d\n", nlive);

  double efr;
  fscanf(input_file, "%lf", &efr);
  printf("efr is %lf\n", efr);

  double tol;
  fscanf(input_file, "%lf", &tol);
  printf("tol is %lf\n", tol);

  int ndims;
  fscanf(input_file, "%d", &ndims);
  printf("ndims is %d\n", ndims);

  int nPar;
  fscanf(input_file, "%d", &nPar);
  printf("nPar is %d\n", nPar);

  int nClsPar;
  fscanf(input_file, "%d", &nClsPar);
  printf("nClsPar is %d\n", nClsPar);

  int updInt;
  fscanf(input_file, "%d", &updInt);
  printf("updInt is %d\n", updInt);

  double Ztol;
  fscanf(input_file, "%lf", &Ztol);
  printf("Ztol is %lf\n", Ztol);

  int maxModes;
  fscanf(input_file, "%d", &maxModes);
  printf("maxModes %d\n", maxModes);

  int seed;
  fscanf(input_file, "%d", &seed);
  printf("seed is %d\n", seed);

  int resume;
  fscanf(input_file, "%d", &resume);
  printf("resume is %d\n", resume);

  int outfile;
  fscanf(input_file, "%d", &outfile);
  printf("outfile is %d\n", outfile);

  int maxiter;
  fscanf(input_file, "%d", &maxiter);
  printf("maxiter is %d\n", maxiter);

  int initMPI;
  fscanf(input_file, "%d", &initMPI);
  printf("initMPI is %d\n", initMPI);

  fscanf(input_file, "%d", &USE_DD);
  printf("USE_DD is %d\n", USE_DD);
  
  fscanf(input_file, "%d", &USE_ICECUBE);
  printf("USE_ICECUBE is %d\n", USE_ICECUBE);
  printf("\n");
  
  
  //Read in some information from params.ini  
  USE_SI = read_param_int("params.ini", "USE_SI");
  printf("USE_SI is %d\n", USE_SI);

  USE_SD = read_param_int("params.ini", "USE_SD");
  printf("USE_SD is %d\n", USE_SD);

  vmode = read_param_int("params.ini", "vmode");
  printf("vmode is %d\n", vmode);

  int decay_channel = read_param_int("params.ini", "decay_channel");
  printf("decay channel is %d\n", decay_channel);

  //Read in some speed distribution information from dist.txt
  int N_dist = read_param_int("dist.txt", "N_dist");
  double rho = read_param_double("dist.txt", "rho_x");
  
  N_expt = 0;
  FLOAT_BG = 0;

  //if using direct detection, read in other information from params.ini
  if (USE_DD)
    {  
      N_expt = read_param_int("params.ini", "N_expt");        
      FLOAT_BG = read_param_int("params.ini", "USE_FLOAT_BG");  
      USE_VARY_FF = read_param_int("params.ini", "USE_VARY_FF");

      printf("Direct detection parameters:\n");
      printf("\tN_expt is %d\n", N_expt);
      printf("\tFLOAT_BG is %d\n", FLOAT_BG);
      printf("\tUSE_VARY_FF is %d\n", USE_VARY_FF);

    }

  int nSimplexPar;
  if (vmode == 1)
    {
      nSimplexPar =   nPar - USE_SI - USE_SD - 1 - N_expt*FLOAT_BG;
    }
  else
    {
      nSimplexPar = 0;
    }

  //-----------------------------------
  nSimplexPar = 0;
 
  printf("nSimplexPar is %d\n", nSimplexPar);
  printf("\n");

  ScalesMin = (double*)malloc(ndims * sizeof(double));
  ScalesMax = (double*)malloc(ndims * sizeof(double));
  for (i = 0; i < ndims; i++)
    {
      fscanf(input_file, "%lf", &ScalesMin[i]);
      fscanf(input_file, "%lf", &ScalesMax[i]);
      printf("Range for parameter %d is %lf %lf\n", i, ScalesMin[i], ScalesMax[i]);
    }
  printf("\n");


  int pWrap[ndims];
  for(i = 0; i < ndims; i++)
    pWrap[i] = 0;

  fclose(input_file);

  int unphys = 0;
  int hwarning = 0;
  int acceptable = 0;
  int iend = 0;
  int ierr = 0;

  dsinterface_init_(&rho, &N_dist, &decay_channel);

  printf("Ready to run\n");

  printf("NOTE: Currently using reweighting of bins...\n");

  printf("NOTE: Floating nuclear physics is not correctly supported...\n");

  FILE *chain_file = fopen(chain_file_name,"r");

  FILE *chain_file_out = fopen(chain_file_name_output, "w");

  double item;
  char status;

  double marg;
  double lk;

  double* params;

  params = (double*)malloc(ndims*sizeof(double));

  int counter = 0;

  double* Ne;
  Ne = (double*)malloc((N_expt+1)*sizeof(double));

  //Count number of lines in the file
  int N_lines = 0;
  char ch;
  while(!feof(chain_file))
  {
    ch = fgetc(chain_file);
    if(ch == '\n')
    {
      N_lines++;
    }
  }

  fseek(chain_file, 0, SEEK_SET); 

  int line = 0;
  int pcounter = 0;

  while(fscanf(chain_file,"%lf",&item) == 1)  
  {    
       
       //printf("\n%e", item);

       if (counter == 0) marg = item;
       if (counter == 1) lk = item;

       if (counter > 1)
       {
           params[counter-2] = item;
       }
       
       counter++;
       if (counter == (ndims + 2))
       {
           //Do the calculation, print the result and reset the counter...
           counter = 0;
           line++;
           LogLike(params, &ndims, &nPar, Ne);

           fprintf(chain_file_out,"%e\t%e\t", marg, lk);
           int ic;
           for (ic = 0; ic < ndims; ic++)
           {
              fprintf(chain_file_out,"%e\t", params[ic]);
           }
           //printf("\n");
           for (ic = 0; ic < N_expt; ic++)
           {
              fprintf(chain_file_out,"%e\t", Ne[ic]);
           }
           fprintf(chain_file_out, "%e\n", Ne[N_expt]);
           //fprintf(chain_file_out,"\n");

           double frac = 100*line/N_lines;
           if (frac > pcounter)
           {
              printf("%d percent completed...\n", pcounter);
              pcounter++;
           }
       }
  }    


  fclose(chain_file);
  fclose(chain_file_out);

  //run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, nSimplexPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLike, dumper, context);

  free(Ne);
  free(params);

  free(ScalesMax);
  free(ScalesMin);

  return 0;
}


