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

// Project headers
#include "../wham/GptlWrappers.h"
#include "../wham/WhamDriver.h"

int main(int argc, char* argv[]) 
{
	// Input checking
	if ( argc < 2 ) {
		std::cerr << "WHAM: Error - missing input file";
		return 1;
	}

	std::string options_file(argv[1]);

	//----- Run Wham -----//

	// GPTL: start
	GPTL::GlobalOptions options;  // sets default options
	GPTL::initialize();           // Initialize GPTL
	GPTL::Timer timer("WHAM");
	timer.start();

	WhamDriver wham_driver(options_file);
	wham_driver.run_driver();

	// GPTL: done
	timer.stop();
	GPTL::print("gptl.log");

	return 0;
}
