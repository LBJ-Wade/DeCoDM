#ifndef EVENTCLASS_H
#define EVENTCLASS_H

#include <fstream>

class Event
{
  public:
    
    double energy; //recoil energy in keV
    double theta;  //in radians
    double phi;    //in radians
    double time;   //in seconds from 2000 (MJD 51544.000000000000)
    
    
    
    //Constructors
    Event(double E)
	    :energy(E)
    {
        theta = 0;
        phi = 0;
        time = 0;
    };
    
    Event(double E,double t,double p)
    :energy(E), theta(t), phi(p)
    {
        time = 0;
    };
    Event(double E,double t,double p, double ti)
    :energy(E), theta(t), phi(p), time(ti){};
    
    void displayEvent(); //Print to screen
    void printEvent(std::ofstream* outputfile); //Print Event onto a line in outputfile
    
};


#endif
