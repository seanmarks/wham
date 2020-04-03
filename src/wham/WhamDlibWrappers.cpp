#include "WhamDlibWrappers.h"

//----- Objective Function Wrapper -----//

WhamDlibEvalWrapper::WhamDlibEvalWrapper(Wham& wham) 
 : wham_ptr_(&wham) 
{
};


double WhamDlibEvalWrapper::operator() (const Wham::ColumnVector& df) const
{
	return wham_ptr_->evalObjectiveFunction(df);
};



//----- Objective Function Derivatives Wrapper -----//

WhamDlibDerivWrapper::WhamDlibDerivWrapper(Wham& wham) 
 : wham_ptr_(&wham) 
{
}


const Wham::ColumnVector WhamDlibDerivWrapper::operator() (const Wham::ColumnVector& df) const
{
	return wham_ptr_->evalObjectiveDerivatives(df);
};
