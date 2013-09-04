#ifndef EVENTCLASS_H
#define EVENTCLASS_H

#include <fstream>

class Event
{
  public:
    
    double energy; //recoil energy in keV
    double theta;  //in radians
    double phi;    //in radians
    
    //Constructor
    Event(double E,double t,double p)
	    :energy(E), theta(t), phi(p) {};
    
    
    void displayEvent(); //Print to screen
    void printEvent(std::ofstream* outputfile); //Print Event onto a line in outputfile
    
};


#endif
