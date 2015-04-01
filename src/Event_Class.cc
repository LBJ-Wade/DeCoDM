#ifndef EVENT_CLASS_CC
#define EVENT_CLASS_CC

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Event_Class.h"

//Print to screen
void Event::displayEvent()
{
  std::cout << "E (keV): "<< energy << "\t theta (rad): " << theta << "\t phi (rad): " << phi << std::endl;
}

//Print Event into a line in outputfile
void Event::printEvent(std::ofstream* outputfile)
{
  *outputfile << energy << "\t" << theta << "\t" << phi << std::endl;
}

#endif
