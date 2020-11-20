// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "OrderParameter.h"

//#include "Assert.hpp"

OrderParameter::OrderParameter(
	const std::string& name,
	const ParameterPack& input_pack
):
	name_(name),
	bins_()
{
	using KeyType = ParameterPack::KeyType;

	// Histogram settings
	const ParameterPack* bins_pack_ptr = input_pack.findParameterPack("Bins", KeyType::Required);
	bins_.setBins( *bins_pack_ptr );
}