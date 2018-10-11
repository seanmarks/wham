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
	if ( argc < 5 ) {
		std::cerr << "Usage:  ./WHAM.exe <options_file> <data_summary> "
		          << "<biasing_params_file> <time_series_y_files_list>";
		return 1;
	}

	std::string options_file(argv[1]);
	std::string data_summary_file(argv[2]);
	std::string biasing_parameters_file(argv[3]);
	std::string time_series_y_files_list(argv[4]);

	//----- Run Wham -----//

	Wham wham(options_file, data_summary_file, biasing_parameters_file,
	          time_series_y_files_list);
	wham.printRawDistributions();
	wham.solve();

	return 0;
}
