/*	main.cpp
 *
 *	ABOUT:
 *	AUTHOR: Sean M. Marks
 *	UPDATED: 
 *	NOTES:
 *	DEVELOPMENT:
 */

// Standard headers
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Library headers
#include "dlib/optimization.h"

// Project headers
#include "Wham.h"

int main(int argc, char* argv[]) 
{
	// Input checking
	if ( argc < 2 ) {
		std::cerr << "Usage:  ./WHAM.exe <options_file>";
		return 1;
	}

	std::string options_file(argv[1]);

	//----- Run Wham -----//

	Wham wham(options_file);
	wham.run_driver();

	return 0;
}
