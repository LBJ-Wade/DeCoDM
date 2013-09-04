#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>
#include <string.h>
//#include <gsl/gsl_interp.h>
#include "RadonTransform.h"


#define PI 3.14159265358
#define deltaV 25
#define deltaQ 25


//At the moment everything is in degrees!!!
//i,j are zero based indices
//Currently, q is an 'integer' index...change so that I specify at different values...

//fout isn't necessarily going to have the same number of points as fin...

//Is there a VERY discretised matrix method I can use to evaluate the radon transform

//Filtered back projection

//May need a window function...

//Pixel inversion...

//Event generation

//Event rate calculation...

//Maybe do some interpolation so that I can specify Nq...

//--------Function Prototypes----------

void rotateImage(double* Rout, double* Rin, double theta, int Nx, int Ny);
double rotatePointX(int i, int j, double theta);
double rotatePointY(int i, int j, double theta);
//int getIndex(int i, int j, int Nx);
double pointRadonTransform(double* fin, int Nx, int Ny, double q, double theta);
void imageRadonTransform(double* fout, double* fin, int Nx, int Ny, int Ntheta);
void imageRadonTransform(double* fout, double* fin, int Nx, int Ny);
void maxwellDistribution(double* f, int Nx, int Ny, double vlagX, double vlagY, double sigma);
void squareDistribution(double* f, int Nx, int Ny, double vlagX, double vlagY, int size);
void inverseRadonTransform(double* fout, double* fin, int Nx, int Ny, int Ntheta);
void applyRampFilter(double* fout, double* fin, int Nx, int Ny);
//int getIndexRowMajor(int i, int j, int Ny);
void addNoise(double* f, int Nx, int Ny, double noise);
void addPointNoise(double* f, int Nx, int Ny, double probability, double size);
//void printMatrix(double* f, char* filename, int Nx, int Ny, int Nz);


void squareDistribution(double* f, int Nx, int Ny, int Nz, double vlagX, double vlagY, double vlagZ, int size);
void addNoise(double* f, int Nx, int Ny, int Nz, double noise);
void addPointNoise(double* f, int Nx, int Ny, int Nz, double probability, double size);
void mask(double *f, int Nq, int Ntheta, int Nphi, double qmin, double qmax);
void realmask(double *f, int Nx, int Ny, int Nz, double rmin, double rmax);
void nullDistribution(double* f, int Nx, int Ny, int Nz);
double diffRate (double qx, double qy, double m_n, double m_x);
void generateEvents(double* f, double* rate, void * params, int Nq, int Ntheta, int Nphi, double m_n, double m_x);
double N_expected(double* rate, int Nq, int Ntheta, int Nphi, double m_n, double m_x);
void generateRate(double* rate, int Nq, int Ntheta, int Nphi, double m_n, double m_x);

//--------Global Variables-------------
//gsl_rng * r;


//--------Function Declarations--------

//CHECK THE 2-D STUFF WORKS!

void rotateImage(double* Rout, double* Rin, double theta, int Nx, int Ny)
{

  //Here's we're using an integer representation (i.e. multiply everything by deltaX, deltaY to get proper distances

  //May need to zero-pad Rout

  //May get aliasing problems...

  double *x1;
  double *x2;

  x1 = (double*)malloc(Nx*Ny*sizeof(double));
  x2 = (double*)malloc(Nx*Ny*sizeof(double));

  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);

  //Now loop over and get the rotated coordinates NB: rotate backwards to get coordinates of old image based on new one
  int m = 0;

  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     x1[m] = round(xmid+rotatePointX(i-xmid,j-ymid,-theta));	//Don't round these off - interpolate...
     x2[m] = round(ymid+rotatePointY(i-xmid,j-ymid,-theta));

     //Check that the new points lie in the image and fix if not
     if ((x1[m] < 0)||(x1[m] > Nx)) x1[m] = 0;
     if ((x2[m] < 0)||(x2[m] > Ny)) x2[m] = 0;

     m++;
   }
  }

  //Here - implement gsl bilinear interpolation...

  //Interpolate...do i need linear or bilinear. Only interpolating along q? Maybe just linear for now???

  //Also need GSL for random numbers...

  int n = 0;
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {

      Rout[getIndex(i,j,1,Nx,1)] = Rin[getIndex(x1[n],x2[n],1, Nx,1)];
      //output += (1+x1 - q/deltaQ)*frot[getIndex(i, x1+ymid, Nx)] + (q/deltaQ - x1)*frot[getIndex(i, x2+ymid, Nx)];
      n++;
    }
  }

  free(x1);
  free(x2);

}

//WAIT! - NEED TO TAKE THEM TO XMID!!!

double rotatePointX(int i, int j, double theta)
{
  //Convert theta to radians
  theta *= PI/(180.0);

  //Evaluate new x-coord
  double newX = (i*cos(theta) - j*sin(theta));

  return newX;
}

double rotatePointY(int i, int j, double theta)
{
  //Convert theta to radians
  theta *= PI/(180.0);

  //Evaluate new y-coord
  double newY = (i*sin(theta) + j*cos(theta));

  return newY;
}

int getIndex(int i, int j, int Nx)
{
 return (i + j*Nx);
}

int getIndexRowMajor(int i, int j, int Ny)
{
 return (j + i*Ny);
}

double pointRadonTransform(double* fin, int Nx, int Ny, double q, double theta)
{
 double output = 0;

 int ymid = ceil((Ny-1.0)/2.0);

 double *frot;
 frot= (double*)malloc(Nx*Ny*sizeof(double));
 rotateImage(frot, fin, theta, Nx, Ny);

 //Initialise interpolation objects
  //gsl_interp * interp = gsl_interp_alloc (gsl_interp_linear, Nx)

  //I've got backwards and forwards acting differently here...MAKE SURE I TREAT THEM CORRECTLY!!!

  int x1 = 0;
  int x2 = 0;

  //if (q >= 0)
  //{
   x1 = floor(q/deltaQ);
   x2 = ceil(q/deltaQ);
  //}
  //else
  //{
   //x1 = ceil(q/deltaQ);
   //x2 = floor(q/deltaQ);
  //}

 for (int i = 0; i < Nx; i++)
 {
  output += (1+x1 - q/deltaQ)*frot[getIndex(i, x1+ymid, Nx)] + (q/deltaQ - x1)*frot[getIndex(i, x2+ymid, Nx)];
  //output += frot[getIndex(i, round(q/deltaQ)+ymid, Nx)];
 }

 free(frot);
 //gsl_interp_free(interp);

 return output;
}
/*
void imageRadonTransform(double* fout, double* fin, int Nx, int Ny, int Ntheta)
{
  //What shape should the new fout be - can we decide beforehand which theta and q values we want. YES!!! For now Nq = Ny...and use
  //fout has shape Ntheta x Nq
  int Nq = Ny;
    int ymid = ceil((Ny-1.0)/2.0);
 double *frot;
 frot = (double*)malloc(Nx*Ny*sizeof(double));

 for (int k = 0; k < Ntheta; k++)
 {
    rotateImage(frot, fin, k*180.0/Ntheta, Nx, Ny);
    for (int q = 0; q < Nq; q++)
    {
      int index = getIndex(k,q,Ntheta);
      fout[index] = 0;
      for (int i = 0; i < Nx; i++)
      {
	fout[index] += frot[getIndex(i, q, Nx)];
      }
    }
 }

 free(frot);
}
*/

void imageRadonTransform(double* fout, double* fin, int Nx, int Ny) //Speed this up, its rather slow...
{									//Allow different Nx, Ny for the two in and out...? Why???
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  double q = 0;
  double theta = 0;


  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     q = (sqrt(pow((i-xmid),2) + pow((j-ymid),2)));
    if (i == xmid)
    {
      theta = 0;
  //    q = (j - ymid);
    }
    else
    {
      theta = 180*atan(1.0*(j-ymid)/(i-xmid))/(PI);
    }
    if ((i  < xmid))
    {
      //q = -q;
    }
    fout[getIndex(i,j,Nx)] = 0;
    if (abs(q) <= Nx/2)
    {
      fout[getIndex(i,j,Nx)] = pointRadonTransform(fin, Nx, Ny, q, theta);
    }
   }
  }
}

void maxwellDistribution(double* f, int Nx, int Ny, double vlagX, double vlagY, double sigma)
{
    //Why is this offset - fix it - zero indexed...

    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);

    double Norm = pow(2*PI*sigma*sigma,-3.0/2.0); ///////////Normalisation for 2D or 3D

    double Norm2 = 0;

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] = Norm*exp(-0.5*(pow(((i - xmid)*deltaV - vlagX),2) + pow((j-ymid)*deltaV - vlagY,2))/(sigma*sigma));
       Norm2 += f[getIndex(i,j,Nx)];
     }
    }

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] /= Norm2;

     }
    }

}



void squareDistribution(double* f, int Nx, int Ny, double vlagX, double vlagY, int size)
{
    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);

    //This method is off by a bit...

    double Norm;

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] = 0;
       if (((i - xmid) >= (vlagX -size/2))&&((i - xmid) <= (vlagX + size/2)))
       {
	        if (((j - ymid) > (vlagY -size/2))&&((j - ymid) < (vlagY + size/2)))
		{
		         f[getIndex(i,j,Nx)] = 1;
			 Norm++;
		}
       }

     }
    }

     for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] /= Norm;

     }
    }

}

void inverseRadonTransform(double* fout, double* fin, int Nx, int Ny, int Ntheta)
{


  //Gondolo method B.4
  //First evaluate the integral at each point...

  //There is a slight offset in the reconstruction...in v-space - maybe due to fudgy rotations....

  //Get this working properly - there is some angular abberation for some reason
  //Then think about applying the RAMP filter...

    //If I'm ramp filtering it might make sense to use the Fourier slice method...

  //Initialise 'intermediate' ramped transform
  double *f;


  f = (double*)malloc(Ntheta*Ny*sizeof(double));

  //std::cout << "Before filtering..." << std::endl;

  applyRampFilter(f,fin,Ntheta,Ny);
  //also change fin down there

  //std::cout << "After filtering..." << std::endl;

  //Possibly need some midpoints...
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);

  double Norm = 0;

  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     fout[getIndex(i,j,Nx)] = 0;
     for (int t = 0; t < Ntheta; t++)
     {
       //Check the integration angles...
       //May need interpolation here also...
       //Possibly a slight problem with directionality...
	double theta =  PI/2.0 + t*PI/Ntheta; //Also rotate by pi/2...for some reason

	//Do I need factors of deltaQ and deltaV here???
	double q = -round(((i-xmid)*cos(theta)+(j-ymid)*sin(theta)));

	       //std::cout << t <<"\t" << theta << "\t" << q << std::endl;
	if ((abs(q) < Ny/2))
	{
	  //Where do these factors come from...?
	  fout[getIndex(i,j,Nx)] += (1.0/(deltaV*deltaV))*f[getIndex(t,q+ymid,Ntheta)]/(1.0*Nx*Ntheta);
	}
     }
     Norm += fout[getIndex(i,j,Nx)];
    }
  }

  //Somewhere I'm out by a funky factor of 2 - I think it's in the original normalisation!
  //Funny inverting fudge...


  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     fout[getIndex(i,j,Nx)] /= Norm;
   }
  }

  free(f);

}

void applyRampFilter(double* fout, double* fin, int Nx, int Ny)
{
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);


  //NO!, I just want to apply the FFT along the q direction, not the theta direction. Does it make a difference? No.

  fftw_complex *in;
  fftw_complex *out;
  fftw_plan pForward;
  fftw_plan pBackward;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
  pForward = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pBackward = fftw_plan_dft_2d(Nx, Ny, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);


  //Pack the input array
   for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      in[getIndexRowMajor(i,j,Ny)][0] = fin[getIndex(i,j,Nx)];
      in[getIndexRowMajor(i,j,Ny)][1] = 0;
    }
  }




  //std::cout << "Packing done..." << std::endl;

  //Execute forward FFT
  fftw_execute(pForward);



  //std::cout << "FFT done..." << std::endl;

  //Apply ramp filtering - fudge so that the ramp is even...NEED TO APPLY THIS PROPERLY...
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < ymid; j++)
    {
      out[getIndexRowMajor(i,j,Ny)][0] *= deltaQ*(j);
      out[getIndexRowMajor(i,j,Ny)][1] *= deltaQ*(j);
      //Apply window function
     // out[getIndexRowMajor(i,j,Ny)][0] *= 0.5*(1+cos(2.0*PI*j/(Ny/2-1.0)));
     // out[getIndexRowMajor(i,j,Ny)][1] *= 0.5*(1+cos(2.0*PI*j/(Ny/2-1.0)));
    }
    for (int j = ymid; j < Ny; j++)
    {
     out[getIndexRowMajor(i,j,Ny)][0] *= deltaQ*abs(Ny-j);
     out[getIndexRowMajor(i,j,Ny)][1] *= deltaQ*abs(Ny-j);
     //Apply window functions
     //out[getIndexRowMajor(i,j,Ny)][0] *= 0.5*(1+cos(2.0*PI*abs(Ny-j)/(Ny/2-1.0)));
     //out[getIndexRowMajor(i,j,Ny)][1] *= 0.5*(1+cos(2.0*PI*abs(Ny-j)/(Ny/2-1.0)));
    }
  }

  //Add the DC offset

  //If I offset the 'j's I get better results - a discretisation error - need interpolation??
  //Offset of about +0.2 works well in general...Aaaah but this is because you end up with negative bits...

  //std::cout << "Filter applied..." << std::endl;

  //Execute backward FFT
  fftw_execute(pBackward);

  //std::cout << "IFFT done..." << std::endl;

  //Unpack FFT
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      fout[getIndex(i,j,Nx)] = in[getIndexRowMajor(i,j,Ny)][0];
    }
  }
  for (int i = 0; i < Nx*Ny; i++)
  {

   //std::cout << in[i][0] << std::endl;
   //std::cout << in[i][1] << std::endl;
   //in[i][1] = 0;
  }

  //std::cout << "Unpacking done..." << std::endl;

  //Clear memory
  fftw_destroy_plan(pForward);
  fftw_destroy_plan(pBackward);
  fftw_free(in); fftw_free(out);
}

//NB: rescaled by number of points in grid...
void addNoise(double* f, int Nx, int Ny, double noise)
{
    gsl_rng * r;
    const gsl_rng_type * T;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,(unsigned)time(NULL));

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
      f[getIndex(i,j,Nx)] += gsl_ran_poisson(r,noise)*1.0/(Nx*Ny);
     }
    }


}

void addPointNoise(double* f, int Nx, int Ny, double probability, double size)
{
    gsl_rng * r;
    const gsl_rng_type * T;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,(unsigned)time(NULL));

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       if (gsl_ran_flat(r,0,1) < probability)
       {
	  f[getIndex(i,j,Nx)] += size;
       }
     }
    }


}

//Can probably improve with interpolation...
/*
void printMatrix(double* f, char* filename, int Nx, int Ny, int Nz) //Extend to 3-D
{
    std::ofstream outputfile;
    outputfile.open (filename);
    for (int j = 0; j < Ny; j++)
    {
      for (int i = 0; i < Nx; i++)
      {
	outputfile << f[getIndex(i,j,Nx)] << "\t";
      }
      outputfile << "\n";
    }
    outputfile.close();

}
*/
//Definitely need to do some interpolation...


//-----------------------------------------------------------------------------
//-----------Left over code from actually implementing these calculations------
//-----------------------------------------------------------------------------


void mask(double *f, int Nq, int Ntheta, int Nphi, double qmin, double qmax)
{
    for (int i = 0; i < Nq; i++)
    {
      double q =  deltaQ*(i-(ceil((Nq-1.0)/2)));
      if ((abs(q) < qmin)||(abs(q) > qmax))
      {
	for (int j = 0; j < Ntheta; j++)
	{
	  for (int k = 0; k < Nphi; k++)
	  {
	    f[getIndex(i,j,k,Nq,Ntheta)] = 0;
	  }
	}
      }
    }


}

void realmask(double *f, int Nx, int Ny, int Nz, double rmin, double rmax)
{
    int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  int zmid = ceil((Nz-1.0)/2.0);

    for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
	{
	  for (int k = 0; k < Nz; k++)
	  {
	    double r =  deltaV*sqrt(pow(i-xmid,2) + pow(j-ymid,2) + pow(k-zmid,2));
	      double theta = acos(deltaV*(k-zmid)/r);
	    if ((r < rmin)||(r > rmax))
	    {

	      f[getIndex(i,j,k,Nx,Ny)] = 0;
	    }
	    if ((r > rmax))
	    {
	      //f[getIndex(i,j,k,Nx,Ny)] -= 2e-4*pow(r,-5.0/2.0)*0.5*(3*cos(theta)*cos(theta)-1);
	    }
	}
      }
    }

}


//Need two forms - a pixel one and an angle one...







void generateRate(double *rate, int Nq, int Ntheta, int Nphi, double m_n, double m_x)
{
  double theta,phi,v;


  for (int i = 0; i < Nq; i++)
  {
    for (int t = 0; t < Ntheta; t++)
    {
      for (int p = 0; p < Nphi; p++)
      {
	    theta = PI/2 + 0.5*t*PI/Ntheta;
	    phi = 2*p*PI/Nphi;
	    v = deltaQ*(i-ceil((Nq-1.0)/2));

	    //rate[getIndex(i,t,p,Nq,Ntheta)] = diffRate(v, theta, phi, m_n, m_x);
      }
    }
  }

}

//-------------------------------MAKE SURE MY LIKELIHOOD IS REASONABLE---------------------------------------
/*
double N_expected(double* rate, int Nq, int Ntheta, int Nphi, double m_n, double m_x)
{
  double theta,v;
  double Ne = 0;



 for (int i = 0; i < Nq; i++)
 {
   v = deltaQ*(i-ceil((Nq-1.0)/2));
   //Exclude the vq = 0
   if (v != 0)
   {
      double deltaE =  2*abs(v_min_inverse(v, m_n, m_x))*deltaQ/abs(v); //Approximately...?
      for (int t = 0; t < Ntheta; t++)
      {
	theta = PI/2 + 0.5*t*PI/Ntheta;
	double deltaTheta = 0.5*PI/Ntheta;
	for (int p = 0; p < Nphi; p++)
	{
	  double deltaPhi = PI/Nphi;
	  if (v_min_inverse(v, m_n, m_x)> 0)	//Threshold!!!
	  {
	      Ne += rate[getIndex(i,t,p,Nq,Ntheta)]*deltaE*sin(theta)*deltaTheta*deltaPhi;
	  }
	}
      }
   }
 }
 return Ne;
}
*/

void generateEvents(double* f, double* rate, void* params, int Nq, int Ntheta, int Nphi,double m_n, double m_x)
{
  //Generate events and store as pixel intensities
  gsl_rng * r;
  const gsl_rng_type * T;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,(unsigned)time(NULL));

    //double Ne = N_expected(maxwell,params);

    std::cout << "Wrong version " << std::endl;

    double Ne = 10;

    int No = gsl_ran_poisson(r,Ne);

    std::cout << "Number expected: " << Ne << "; Number observed: " << No << std::endl;

    for (int N = 0; N < No; )
    {

      //-------------------------------------------------------------------------------------------------------------------
      //-----------This is where the work needs to be done - need to convert between all q and positive q properly... (should be easy)
      //-----------Then need to clean up and document code..>!!!!!!!!!!!!!!
      //-------------------------------------------------------------------------------------------------------------------

      double t = gsl_rng_uniform(r);
      double p = gsl_rng_uniform(r);
      double q = gsl_rng_uniform(r);

      if (gsl_ran_flat(r,0,1) < rate[getIndex(q,t,p,Nq,Ntheta)])
      {
	f[getIndex(q,t,p,Nq,Ntheta)] += 1.0;
	N++;
      }

    }
}



//Need to add some 'standard' halos

void maxwellDistribution(double* f, int Nx, int Ny, int Nz, double vlagX, double vlagY, double vlagZ, double sigma)
{
    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);
    int zmid = ceil((Nz-1.0)/2.0);

    double Norm = pow(2*PI*sigma*sigma,-3.0/2.0);

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       for (int k = 0; k < Nz; k++)
       {
	 double temp;
	 temp = Norm*exp(-0.5*(pow(((i - xmid)*deltaV - vlagX),2) +
						pow(((j - ymid)*deltaV - vlagY),2) +
						pow(((k - zmid)*deltaV - vlagZ),2)
						)/(sigma*sigma));
	  f[getIndex(i,j,k,Nx,Ny)] = temp;
	  //Norm2 += temp;
       }
     }
    }

}

void assDistribution(double* f, int Nx, int Ny, int Nz)
{
    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);
    int zmid = ceil((Nz-1.0)/2.0);

    double sigma = 156;
    double vlag = 230;



    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       for (int k = 0; k < Nz; k++)
       {
	 double v = sqrt(pow(((i - xmid)*deltaV),2) +
		      pow(((j - ymid)*deltaV),2) +
		      pow(((k - zmid)*deltaV),2) + 1);
	//std::cout << v << std::endl;
	 double temp = 0;
	 if (k > zmid)
	 {
	   temp = (0.5/(sqrt(2*PI)*vlag*sigma*v))*exp(-(v*v + vlag*vlag)/(2*sigma*sigma))*(exp(v*vlag/(sigma*sigma)) - 1);
	   //temp = 1;
	 }
	 else
	 {
	   temp = (0.5/(sqrt(2*PI)*vlag*sigma*v))*exp(-(v*v + vlag*vlag)/(2*sigma*sigma))*(1 - exp(-v*vlag/(sigma*sigma)));
	   //temp = 0;
	 }
	 //std::cout << temp << std::endl;
	   f[getIndex(i,j,k,Nx,Ny)] = temp;
	  //Norm2 += temp;
       }
     }
    }

}

void nullDistribution(double* f, int Nx, int Ny, int Nz)
{
    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);
    int zmid = ceil((Nz-1.0)/2.0);

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       for (int k = 0; k < Nz; k++)
       {
	 double x = deltaV*1.0*(i-xmid);
	 double y = deltaV*1.0*(j-ymid);
	 double z = deltaV*1.0*(k-zmid);
;
	 double r = sqrt(x*x + y*y + z*z);
	  double theta = acos(z/r);
	 //std::cout << r << "\t" << theta << std::endl;
	 if (((x!=0)||(y!=0)||(z!=0)))
	 {
	   if (r < 450)
	   {
	     //f[getIndex(i,j,k,Nx,Ny)] +=2*(1.1284 + 1.5045*((1/r*r)-1))*sin(2*theta)/(r*r);
	    f[getIndex(i,j,k,Nx,Ny)] += (1e-5)*0.5*(3*cos(theta)*cos(theta)-1)/(r*r*r);
	    //f[getIndex(i,j,k,Nx,Ny)] += (1e-4)*(35*pow(cos(theta),4) - 30*pow(cos(theta),2) + 3)/(r*r*r);
	    //f[getIndex(i,j,k,Nx,Ny)] += (10000e-9)*(5*pow(cos(theta),3) - 3*cos(theta))/(pow(r,4));
	  }
	 }
	//f[getIndex(i,j,k,Nx,Ny)] += 0.001;
       }
     }
    }

}

void squareDistribution(double* f, int Nx, int Ny, int Nz, double vlagX, double vlagY, double vlagZ, int size)
{
  /*
    int xmid = ceil((Nx-1.0)/2.0);
    int ymid = ceil((Ny-1.0)/2.0);

    //This method is off by a bit...

    double Norm;

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] = 0;
       if (((i - xmid) >= (vlagX -size/2))&&((i - xmid) <= (vlagX + size/2)))
       {
	        if (((j - ymid) > (vlagY -size/2))&&((j - ymid) < (vlagY + size/2)))
		{
		         f[getIndex(i,j,Nx)] = 1;
			 Norm++;
		}
       }

     }
    }

     for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       f[getIndex(i,j,Nx)] /= Norm;

     }
    }
  */
}

void addNoise(double* f, int Nx, int Ny, int Nz, double noise)
{
    gsl_rng * r;
    const gsl_rng_type * T;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,(unsigned)time(NULL));

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       for (int k = 0; k < Nz; k++)
       {
	  f[getIndex(i,j,k,Nx,Ny)] += gsl_ran_poisson(r,noise);
       }
     }
    }


}

void addPointNoise(double* f, int Nx, int Ny, int Nz, double probability, double size)
{
   gsl_rng * r;
    const gsl_rng_type * T;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,(unsigned)time(NULL));

    for (int i = 0; i < Nx; i++)
    {
     for (int j = 0; j < Ny; j++)
     {
       for (int k = 0; k < Nz; k++)
       {
	  if (gsl_ran_flat(r,0,1) < probability)
	  {
	    f[getIndex(i,j,k,Nx,Ny)] += size;
	  }
       }
     }
    }


}


//Working reasonably well - only shows a weak dependence on the number of pixels...

//At the moment everything is in degrees!!!
//i,j are zero based indices
//Currently, q is an 'integer' index...change so that I specify at different values...

//Pixel inversion...

//Event generation

//Event rate calculation...

//Maybe do some interpolation so that I can specify Nq...

//Maybe redo to allow in-place rotations and transforms...

//Use vectors...





//--------Function Declarations--------


void rotateImage(double* Rout, double* Rin, double theta, int Nx, int Ny, int Nz, int axis)
{
  //May need to zero-pad Rout

  //May get aliasing problems...

   for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
    for (int k = 0; k < Nz; k++)
    {
      int index = getIndex(i,j,k,Nx,Ny);
      if (isnan(Rin[index])) std::cout << index << std::endl;
    }
   }
  }


  double* Rtemp;

  Rtemp = (double*)malloc(Nx*Ny*Nz*sizeof(double));

  double *x;
  double *y;

  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  int zmid = ceil((Nz-1.0)/2.0);

  if (axis == 1)
  {

      x = (double*)malloc(Nx*Ny*sizeof(double));
      y = (double*)malloc(Nx*Ny*sizeof(double));

      //Now loop over and get the rotated coordinates NB: rotate backwards to get coordinates of old image based on new one
    int m = 0;

    for (int i = 0; i < Nx; i++)
    {

      for (int j = 0; j < Ny; j++)
      {
	  x[m] = (xmid+rotatePointX(i-xmid,j-ymid,-theta));
	  y[m] = (ymid+rotatePointY(i-xmid,j-ymid,-theta));
	  //Check that the new points lie in the image and fix if not
	  if ((x[m] < 0)||(ceil(x[m]) >= Nx)) x[m] = 0;
	  if ((y[m] < 0)||(ceil(y[m]) >= Ny)) y[m] = 0;

	  m++;
      }
    }

    int n = 0;
    for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
      {
	for (int k = 0; k < Nz; k++)
	{
	    int x2 = ceil(x[n]);
	    int x1 = floor(x[n]);

	    int y2 = ceil(y[n]);
	    int y1 = floor(y[n]);

	    if (x1 == x2)
	    {
	       if (y1 == y2)
	      {
		Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(x1,y1,k,Nx,Ny)];
	      }
	      else
	      {
		Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(x1,y1,k,Nx,Ny)]*(y2-y[n]) + Rin[getIndex(x1,y2,k,Nx,Ny)]*(y[n]-y1);
	      }
	    }
	    else if (y1 == y2)
	    {
	      Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(x1,y1,k,Nx,Ny)]*(x2-x[n]) + Rin[getIndex(x2,y1,k,Nx,Ny)]*(x[n]-x1);
	    }
	    else
	    {
	      Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(x1,y1,k,Nx,Ny)]*(x2-x[n])*(y2-y[n]) + Rin[getIndex(x2,y1,k,Nx,Ny)]*(x[n]-x1)*(y2-y[n]) +
						Rin[getIndex(x1,y2,k,Nx,Ny)]*(x2-x[n])*(y[n]-y1) + Rin[getIndex(x2,y2,k,Nx,Ny)]*(x[n]-x1)*(y[n]-y1);
	    }
	}

	n++;
      }

    }
  }
  else if (axis == 2)
  {
      x = (double*)malloc(Ny*Nz*sizeof(double));
      y = (double*)malloc(Ny*Nz*sizeof(double));

    //Now loop over and get the rotated coordinates NB: rotate backwards to get coordinates of old image based on new one
    int m = 0;


    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
	  x[m] = (ymid+rotatePointX(j-ymid,k-zmid,-theta));	//Don't round these off - interpolate...
	  y[m] = (zmid+rotatePointY(j-ymid,k-zmid,-theta));
	  //Check that the new points lie in the image and fix if not
	  //Check this for Nx!=Ny!=Nz
	  if ((x[m] < 0)||(ceil(x[m]) >= Ny))
	  {
	   x[m] = 0;
	  }
	  if ((y[m] < 0)||(ceil(y[m]) >= Nz))
	  {
	    y[m] = 0;

	  }

	  m++;
      }
    }

    int n = 0;
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
	for (int i = 0; i < Nx; i++)
	{
	    int x2 = ceil(x[n]);
	    int x1 = floor(x[n]);

	    int y2 = ceil(y[n]);
	    int y1 = floor(y[n]);

	    if (x1 == x2)
	    {
	       if (y1 == y2)
	      {
		Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(i,x1,y1,Nx,Ny)];
	      }
	      else
	      {
		Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(i,x1,y1,Nx,Ny)]*(y2-y[n]) + Rin[getIndex(i,x1,y2,Nx,Ny)]*(y[n]-y1);
	      }
	    }
	    else if (y1 == y2)
	    {
	      Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(i,x1,y1,Nx,Ny)]*(x2-x[n]) + Rin[getIndex(i,x2,y1,Nx,Ny)]*(x[n]-x1);
	    }
	    else
	    {

	      Rtemp[getIndex(i,j,k,Nx,Ny)] = Rin[getIndex(i,x1,y1,Nx,Ny)]*(x2-x[n])*(y2-y[n]) + Rin[getIndex(i,x2,y1,Nx,Ny)]*(x[n]-x1)*(y2-y[n]) +
						Rin[getIndex(i,x1,y2,Nx,Ny)]*(x2-x[n])*(y[n]-y1) + Rin[getIndex(i,x2,y2,Nx,Ny)]*(x[n]-x1)*(y[n]-y1);
	    }
	}
	n++;
      }

    }
  }

  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
    for (int k = 0; k < Nz; k++)
    {
      int index = getIndex(i,j,k,Nx,Ny);
      Rout[index] = Rtemp[index];
    }
   }
  }


  free(x);
  free(y);
  free(Rtemp);
}




double pointRadonTransform(double* fin, int Nx, int Ny, int Nz, double q, double theta, double phi)
{

 double output = 0;

 int ymid = ceil((Ny-1.0)/2.0);

 double *frot1;
 double *frot2;

 frot1 = (double*)malloc(Nx*Nz*Nz*sizeof(double));
  frot2 = (double*)malloc(Nx*Ny*Nz*sizeof(double));

 rotateImage(frot1, fin, phi, Nx, Ny, Nz,2);
 rotateImage(frot2, frot1, theta, Nx, Ny, Nz, 1);
 //Initialise interpolation objects
  //gsl_interp * interp = gsl_interp_alloc (gsl_interp_linear, Nx)

  //I've got backwards and forwards acting differently here...MAKE SURE I TREAT THEM CORRECTLY!!!

  int x1 = 0;
  int x2 = 0;

  //if (q >= 0)
  //{
   x1 = floor(q/deltaQ);
   x2 = ceil(q/deltaQ);
  //}
  //else
  //{
   //x1 = ceil(q/deltaQ);
   //x2 = floor(q/deltaQ);
  //}

 //---------------------ChecK!!!

 for (int j = 0; j < Ny; j++)
 {
   for (int k = 0; k < Nz; k++)
   {
      output += (1+x1 - q/deltaQ)*frot2[getIndex(x1+ymid, j,k , Nx, Ny)] + (q/deltaQ - x1)*frot2[getIndex(x2+ymid, j, k, Nx, Ny)];
      //output += frot[getIndex(i, round(q/deltaQ)+ymid, Nx)];
   }
 }

 free(frot1);
 free(frot2);

 //gsl_interp_free(interp);

 return output;
}


void imageRadonTransform2D(double* fout, double* fin, int Nx, int Ny, int Ntheta)
{
    int Nz = 1;
    int Nq = Ny;

    //Initialise intermediate matrices
    double *frot1;

    frot1 = (double*)malloc(Nx*Ny*Nz*sizeof(double));

    for (int t = 0; t < Ntheta; t++)
    {
	//Rotate around z
        rotateImage(frot1, fin, 1*t*180.0/Ntheta, Nx, Ny, Nz, 1);

	//Select q-slice
	for (int q = 0; q < Nq; q++)
	{
	  int index = getIndex(q,t,0,Nq,Ntheta);
	  fout[index] = 0;

	  //Integrate over line
	  for (int j = 0; j < Ny; j++)
	  {
	      fout[index] += frot1[getIndex(j, q, 0, Nx, Ny)];
	  }
	}
    }

   //Free memory
    free(frot1);

}


void imageRadonTransform3D(double* fout, double* fin, int Nx, int Ny, int Nz, int Ntheta, int Nphi)
{

    int Nq = Ny;

    //Initialise intermediate matrices
    double *frot1;
    double *frot2;

    frot1 = (double*)malloc(Nx*Ny*Nz*sizeof(double));
    frot2 = (double*)malloc(Nx*Ny*Nz*sizeof(double));

    for (int p = 0; p < Nphi; p++)
    {
	//Rotate around z
        rotateImage(frot1, fin, 2*p*180.0/Nphi, Nx, Ny, Nz, 1);
	for (int t = 0; t < Ntheta; t++)
	{
	  //Rotate around x
	  rotateImage(frot2, frot1, 0.5*t*180.0/Ntheta, Nx, Ny, Nz, 2);

	  //Select q-slice
	  for (int q = 0; q < Nq; q++)
	  {
	    int index = getIndex(q,t,p,Nq,Ntheta);
	    fout[index] = 0;

	    //Integrate over plane
	    for (int j = 0; j < Ny; j++)
	    {
	      for (int k = 0; k < Nz; k++)
	      {
		fout[index] += deltaV*deltaV*frot2[getIndex(j, q, k, Nx, Ny)];
	      }
	    }
	  }
	}
    }

   //Free memory
    free(frot1);
    free(frot2);
}


double mass(double* f, int Nq, int Ntheta, int Nphi, int theta, int phi)
{
  double m = 0;

 for (int q = 0; q < Nq; q++)
 {
   m += f[getIndex(q,theta,phi, Nq,Ntheta)];
 }

  return m;
}

double CoM(double* f, int Nq, int Ntheta, int Nphi, int theta, int phi)
{
  double m = 0;

    int qmid = ceil((Nq-1.0)/2.0);

 for (int q = 0; q < Nq; q++)
 {
   m += deltaQ*(q-qmid)*(f[getIndex(q,theta,phi, Nq,Ntheta)]);
 }

  return m;
}

void imageRadonTransform(double* fout, double* fin, int Nx, int Ny, int Nz) //Speed this up, its rather slow...
{									//Allow different Nx, Ny for the two in and out...? Why???
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  int zmid = ceil((Nz-1.0)/2.0);

  //May need to switch theta and phi...

  double q = 0;
  double theta = 0;
  double phi = 0;

  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
    for (int k = 0; k < Nz; k++)
    {
      q = (sqrt(pow((i-xmid),2) + pow((j-ymid),2) + pow((k-zmid),2)));
      if (i == xmid)
      {
	theta = 0;
    //    q = (j - ymid);
      }
      else
      {
	phi = 180*atan(1.0*(j-ymid)/(i-xmid))/(PI);
	theta = 180*acos(1.0*(k-zmid)/q)/PI;
      }
      if ((i  < xmid))
      {
	//q = -q;
      }
      fout[getIndex(i,j,k,Nx,Ny)] = 0;
      if (abs(q) <= Nx/2)
      {
	fout[getIndex(i,j,k,Nx,Ny)] = pointRadonTransform(fin, Nx, Ny, Nz, q, theta, phi);
      }
    }
   }
  }
}

//Be aware of aliasing...

void inverseRadonTransform2D(double* fout, double* fin, int Nx, int Ny, int Ntheta)
{
  //Set Nz = 1 for 2D
  int Nz = 1;

  //Initialise 'intermediate' filtered transform
  double *f;


  //Calculate some mid-points
  int Nq = Ny;
  int qmid = ceil((Nq-1.0)/2.0);
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);

  //Apply |w| filter
  f = (double*)malloc(Nx*Ny*sizeof(double));
  //applyFilter2D(f,fin,Nq,Ntheta);


  //Intregrate over a (hemi-)sphere of radius r/2 passing through the origin (r = sqrt(x^2 + y^2 + z^2)
  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     for (int k = 0; k < Nz; k++)
     {
      fout[getIndex(i,j,k,Nx,Ny)] = 0;
      for (int t = 0; t < Ntheta; t++)
      {
	    double theta =  PI/2 + 1*t*PI/Ntheta;

	   //Calculate q
	    double q = (deltaV/deltaQ)*(qmid+((i-xmid)*cos(theta)+(j-ymid)*sin(theta)));


	    int q1 = floor(q);
	    int q2 = ceil(q);


	    //Interpolate values of radon transform along q
	    if ((q1 >= 0)&&(q2 < Nq))
	    {
	      if (q1 == q2)
	      {
		f[getIndex(i,j,k,Nx,Ny)] += fin[getIndex(q1,t,0,Nq,Ntheta)];
	      }
	      else
	      {
		f[getIndex(i,j,k,Nx,Ny)] += (fin[getIndex(q1,t,0,Nq,Ntheta)]*(q2-q) + fin[getIndex(q2,t,0,Nq,Ntheta)]*(q-q1));
	      }
	    }
	}
//	if (fout[getIndex(i,j,k,Nx,Ny)] < 0) fout[getIndex(i,j,k,Nx,Ny)] = 0;
      }
     }
   }
  //fout[getIndex(xmid,ymid,0,Nx,Ny)] *=0.5;

  applyFilter2D(fout, f, Nx, Ny);

   for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     for (int k = 0; k < Nz; k++)
     {
      if (fout[getIndex(i,j,k,Nx,Ny)] < 0) fout[getIndex(i,j,k,Nx,Ny)] = 0;
     }
   }
  }


  //Unit normalise the result
  normaliseMatrix(fout, Nx, Ny, Nz);

  free(f);
}

void inverseRadonTransform3D(double* fout, double* fin, int Nx, int Ny, int Nz, int Ntheta, int Nphi)
{
  //Initialise 'intermediate' filtered transform
  double *f;

  //Calculate some mid-points
  int Nq = Ny;
  int qmid = ceil((Ny-1.0)/2);
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  int zmid = ceil((Nz-1.0)/2.0);

  //Sampling errors and the Nyquist sampling theorem - how big is the problem...?

  //Apply w^2 filter
  f = (double*)malloc(Nq*Ntheta*Nphi*sizeof(double));
  applyFilter3D(f,fin,Nq,Ntheta,Nphi);
  //Laplacian(f,fin,Nq,Ntheta,Nphi);

  //Intregrate over a (hemi-)sphere of radius r/2 passing through the origin (r = sqrt(x^2 + y^2 + z^2)
  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     for (int k = 0; k < Nz; k++)
     {
      fout[getIndex(i,j,k,Nx,Ny)] = 0;
      for (int t = 0; t < Ntheta; t++)
      {
	for (int p = 0; p < Nphi; p++)
	{
	    double theta =  PI/2 + 0.5*t*PI/Ntheta;
	    double phi = PI/2 -2*p*PI/Nphi;

	   //Calculate q
	    double q = ((deltaV/deltaQ)*(qmid+((i-xmid)*cos(phi)*sin(theta)+(j-ymid)*sin(phi)*sin(theta)+(k-zmid)*cos(theta))));

	    int q1 = floor(q);
	    int q2 = ceil(q);


	    //Interpolate values of radon transform along q
	    if ((q1 >= 0)&&(q2 < Nq))
	    {
	      if (q1 == q2)
	      {
		fout[getIndex(i,j,k,Nx,Ny)] += sin(theta)*f[getIndex(q1,t,p,Nq,Ntheta)];
	      }
	      else
	      {
		fout[getIndex(i,j,k,Nx,Ny)] += sin(theta)*(f[getIndex(q1,t,p,Nq,Ntheta)]*(q2-q) + f[getIndex(q2,t,p,Nq,Ntheta)]*(q-q1));
	      }
	    }
	  }
	}
	if (fout[getIndex(i,j,k,Nx,Ny)] < 0) fout[getIndex(i,j,k,Nx,Ny)] = 0;
     }
    }
  }


  //Unit normalise the result
  normaliseMatrix(fout, Nx, Ny, Nz);

  free(f);

}

void Laplacian(double* fout, double* fin, int Nx, int Ny, int Nz)
{
  //Initialise 'intermediate' filtered transform
  //applyFilter(f,fin,Nq,Ntheta,Nphi);



for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     for (int k = 0; k < Nz; k++)
     {
	fout[getIndex(i,j,k,Nx,Ny)] = 0;
     }
   }
  }
double Norm = 0 ;

  for (int i = 1; i < Nx-1; i++)
  {
   for (int j = 1; j < Ny-1; j++)
   {
     for (int k = 1; k < Nz-1; k++)
     {
      fout[getIndex(i,j,k,Nx,Ny)] = (1.0/(deltaV*deltaV))*(fin[getIndex(i+1, j, k,Nx, Ny)] - 2*fin[getIndex(i,j,k,Nx,Ny)] + fin[getIndex(i-1, j, k,Nx, Ny)]);
      fout[getIndex(i,j,k,Nx,Ny)] += (1.0/(deltaV*deltaV))*(fin[getIndex(i, j+1, k,Nx, Ny)] - 2*fin[getIndex(i,j,k,Nx,Ny)] + fin[getIndex(i, j-1, k,Nx, Ny)]);
      fout[getIndex(i,j,k,Nx,Ny)] += (1.0/(deltaV*deltaV))*(fin[getIndex(i, j, k+1,Nx, Ny)] - 2*fin[getIndex(i,j,k,Nx,Ny)] + fin[getIndex(i, j, k-1,Nx, Ny)]);
      Norm += fout[getIndex(i,j,k,Nx,Ny)];
     }
    }
  }

  for (int i = 0; i < Nx; i++)
  {
   for (int j = 0; j < Ny; j++)
   {
     for (int k = 0; k < Nz; k++)
     {
	fout[getIndex(i,j,k,Nx,Ny)] /= (Norm);
     }
   }
  }



}

void applyFilter2D(double* fout, double* fin, int Nx, int Ny)
{
  //Set Nz = 1;
  int Nz = 1;

  //Calculate mid-cells
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);



  //Initialise FFTW variables and parameters
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan pForward;
  fftw_plan pBackward;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  pForward = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pBackward = fftw_plan_dft_2d(Nx, Ny, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);


  //Pack the input array
   for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
	in[getIndexRowMajor(i,j,k,Ny,Nz)][0] = fin[getIndex(i,j,k,Nx,Ny)];
	in[getIndexRowMajor(i,j,k,Ny,Nz)][1] = 0;
      }
    }
  }

  //Execute forward FFT
  fftw_execute(pForward);

  //Apply ramp filtering (|w|)
  for (int j = 0; j < Ny; j++)
  {
    int yi = j;
    if (j >= ymid) yi = Ny-j;

      for (int k = 0; k < Nz; k++)
      {
	for (int i = 0; i < Nx; i++)
	{
	  int xi = i;
	  if (i >= xmid) xi = Nx - i;

	  out[getIndexRowMajor(i,j,k,Ny,Nz)][0] *= deltaQ*abs(xi);
	  out[getIndexRowMajor(i,j,k,Ny,Nz)][1] *= deltaQ*abs(xi);
	  //Apply window function
	  //out[getIndexRowMajor(i,j,k,Ny,Nz)][0] *= 0.5*(1+cos(2.0*PI*i/(Nx/2-1.0)));
	  //out[getIndexRowMajor(i,j,k,Ny,Nz)][1] *= 0.5*(1+cos(2.0*PI*i/(Nx/2-1.0)));
	}
      }
  }


  //Execute backward FFT
  fftw_execute(pBackward);

  //Unpack FFT
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
      fout[getIndex(i,j,k,Nx,Ny)] = in[getIndexRowMajor(i,j,k,Ny,Nz)][0];
      }
    }
  }

  //Clear memory
  fftw_destroy_plan(pForward);
  fftw_destroy_plan(pBackward);
  fftw_free(in); fftw_free(out);
}


void applyFilter3D(double* fout, double* fin, int Nx, int Ny, int Nz)
{

  //Calculate mid-cells
  int xmid = ceil((Nx-1.0)/2.0);
  int ymid = ceil((Ny-1.0)/2.0);
  int zmid = ceil((Nz-1.0)/2.0);


  //Initialise FFTW variables and parameters
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan pForward;
  fftw_plan pBackward;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);
  pForward = fftw_plan_dft_3d(Nx, Ny, Nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pBackward = fftw_plan_dft_3d(Nx, Ny, Nz, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);


  //Pack the input array
   for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
	in[getIndexRowMajor(i,j,k,Ny,Nz)][0] = fin[getIndex(i,j,k,Nx,Ny)];
	in[getIndexRowMajor(i,j,k,Ny,Nz)][1] = 0;
      }
    }
  }

  //Execute forward FFT
  fftw_execute(pForward);

  //Apply w^2 filtering
for (int j = 0; j < Ny; j++)
  {
    int yi = j;
    if (j >= ymid)
    {
      yi = Ny-j;
    }
      for (int k = 0; k < Nz; k++)
      {
	int zi = k;
	if (k >= zmid)
	{
	  zi = Nz-k;
	}

	for (int i = 0; i < Nx; i++)
	{
	  int xi = i;
	  if (i >= xmid)
	  {
	    xi = Nx - i;
	  }
	  out[getIndexRowMajor(i,j,k,Ny,Nz)][0] *= deltaQ*deltaQ*(xi*xi);
	  out[getIndexRowMajor(i,j,k,Ny,Nz)][1] *= deltaQ*deltaQ*(xi*xi);
	  //Apply window function
	  //out[getIndexRowMajor(i,j,k,Ny,Nz)][0] *= 0.5*(1+cos(2.0*PI*i/(Nx/2-1.0)));
	  //out[getIndexRowMajor(i,j,k,Ny,Nz)][1] *= 0.5*(1+cos(2.0*PI*i/(Nx/2-1.0)));
	}
      }
  }


  //Execute backward FFT
  fftw_execute(pBackward);

  //Unpack FFT
  for (int i = 0; i < Nx; i++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int k = 0; k < Nz; k++)
      {
      fout[getIndex(i,j,k,Nx,Ny)] = in[getIndexRowMajor(i,j,k,Ny,Nz)][0];
      }
    }
  }

  //Clear memory
  fftw_destroy_plan(pForward);
  fftw_destroy_plan(pBackward);
  fftw_free(in); fftw_free(out);
}


