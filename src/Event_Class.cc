#ifndef EVENT_CLASS_CC
#define EVENT_CLASS_CC

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "DMUtils.h"
#include "Event_Class.h"

//Print to screen
void Event::displayEvent()
{	
  std::cout << "E (keV): "<< energy << "\t theta (rad): " << theta << "\t phi (rad): " << phi << std::endl;
  //std::cout << "E_x (keV): "<< E_x << "\t E_y (keV): " << E_y << "\t E_z (keV): " << E_z << std::endl;
}

//Print Event into a line in outputfile
void Event::printEvent(std::ofstream* outputfile)
{
  *outputfile << energy << "\t" << theta << "\t" << phi << std::endl;
}

//Convert from angles to vector components
void Event::Angles_to_Vector()
{
    E_z = energy*sin(theta)*cos(phi);
	E_x = energy*sin(theta)*sin(phi);
	E_y = -energy*cos(theta);
}

//Convert from vector component to angles
void Event::Vector_to_Angles()
{
	energy = sqrt(E_x*E_x + E_y*E_y + E_z*E_z);
	theta = acos(-E_y/energy);
	phi = atan2(E_x, E_z);
	if (phi < 0) phi += 2*PI;
}

#endif
