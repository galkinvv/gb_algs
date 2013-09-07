/**
\file
Реализация не-inline методов CMonomialBase
*/

#include <cassert>
#include "cmonomial.h"
#include "memorymanager.h"

namespace F4MPI{
template<> size_t PODvecSize<CInternalMonomial>::iterator_step=0;

CMonomialBase::Order CMonomialBase::order;
int CMonomialBase::orderParam = 0;
int CMonomialBase::theNumberOfVariables = 0;
int CMonomialBase::degreessize = 0;

CMonomialBase::Order CMonomialBase::getOrder(){
	return order;
}
void CMonomialBase::setOrder(CMonomialBase::Order ord, int anOrderParam){
	order = ord;
	orderParam = anOrderParam;
}

std::string CMonomialBase::varName(int varIdxFrom1, ParserVarNames* names)
{
	if (names)
	{
		int idxInNames = varIdxFrom1-1;
		assert(idxInNames >= 0 && idxInNames < names->size());
		return (*names)[varIdxFrom1-1];
	}
	else
	{
		std::ostringstream res;
		res<<"_X"<<varIdxFrom1;
		return res.str();
	}
}

const std::string getMonomialOrderName(CMonomialBase::Order orderID){
	if (orderID==CMonomialBase::lexOrder) return "lex";
	if (orderID==CMonomialBase::deglexOrder) return "deglex";
	if (orderID==CMonomialBase::degrevlexOrder) return "degrevlex";
	if (orderID==CMonomialBase::blklexOrder) return "blklex";
	return "UNKNOWN";
}

CMonomialBase::Order CMonomialForLabel::orderForLabel = CMonomialBase::degrevlexOrder;
} //namespace F4MPI
