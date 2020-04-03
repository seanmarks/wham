
#pragma once
#ifndef WHAM_DLIB_WRAPPERS_H
#define WHAM_DLIB_WRAPPERS_H

// Standard headers
#include <vector>

// Library headers
#include "dlib/optimization.h"

// Project headers
#include "Wham.h"

class Wham;

//----- Wrapper for WHAM objective function -----//

class WhamDlibEvalWrapper
{
 public:
	using ColumnVector = dlib::matrix<double,0,1>;

	WhamDlibEvalWrapper(Wham& wham);

	double operator() (const ColumnVector& df) const;

 private:
	Wham* wham_ptr_;

};


//----- Wrapper for derivatives of the WHAM objective function -----//

class WhamDlibDerivWrapper
{
 public:
	using ColumnVector = dlib::matrix<double,0,1>;

	WhamDlibDerivWrapper(Wham& wham);

	const ColumnVector operator() (const ColumnVector& df) const;

 private:
	Wham* wham_ptr_;
};

#endif // WHAM_DLIB_WRAPPERS_H
