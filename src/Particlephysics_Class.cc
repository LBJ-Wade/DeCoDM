#include <iostream>
#include <Particlephysics_class.h>

void Particlephysics::PrintAll()
{
	std::cout << "-----Particle physics parameters------" << std::endl;
	std::cout << "\tm_x = " << m_x << " GeV" << std::endl;
	std::cout << "\tsigma_SI = " << sigma_SI << " cm^2" << std::endl;
	std::cout << "\tsigma_SD = " << sigma_SD << " cm^2" << std::endl;
	std::cout << std::endl;
	std::cout << "\tsigma_05 = " << sigma_O5 << " cm^2" << std::endl;
	std::cout << "\tsigma_07 = " << sigma_O7 << " cm^2" << std::endl;
	std::cout << "\tsigma_015 = " << sigma_O15 << " cm^2" << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	
}





