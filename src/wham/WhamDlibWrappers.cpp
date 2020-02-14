#include "WhamDlibWrappers.h"

//----- Objective Function Wrapper -----//

WhamDlibEvalWrapper::WhamDlibEvalWrapper(Wham& wham) 
 : wham_(wham) 
{
};


double WhamDlibEvalWrapper::operator() (const Wham::ColumnVector& df) const
{
	return wham_.evalObjectiveFunction(df);
};



//----- Objective Function Derivatives Wrapper -----//

WhamDlibDerivWrapper::WhamDlibDerivWrapper(Wham& wham) 
 : wham_(wham) 
{
}


const Wham::ColumnVector WhamDlibDerivWrapper::operator() (const Wham::ColumnVector& df) const
{
	return wham_.evalObjectiveDerivatives(df);
};
