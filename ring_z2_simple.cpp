#include "ring_z2_simple.h"
#include "simple_mon.h"
#include "sparse_matrix.h"
#include <algorithm>
#include <deque>
#include <cassert>
#include <unordered_set>
#include <list>

struct SlowMon : SimpleMon {};
struct SlowPolynomial : std::vector<SlowMon> {};

struct RingZ2SimpleBase::Impl {
	//this struct shoulld contain data corresponding to polynomial ring, but not to some polynomials of this ring
	Impl(int keeped_vars_count)
		:keeped_vars_count_(keeped_vars_count)
	{}
	const int keeped_vars_count_;
	std::vector<SlowMon> new_variables;
};

//input polynomials in ring with original monomial count
struct RingZ2SimpleBase::InPolysSetWithOrigMetadata : std::vector<SlowPolynomial> {};

//output polynomials in ring with new monomial count
struct RingZ2SimpleBase::OutPolysSetForVariyingMetadata : std::vector<SlowPolynomial> {};

RingZ2SimpleBase::RingZ2SimpleBase(int var_count):
	impl_(var_count)
{}

bool RingZ2SimpleBase::ReconstructAndInsertNormalizedImpl(const InPolysSetWithOrigMetadata& reconstruction_basis,
        Enumerator<CrossRingInfo::PerVariableData> top_info,
        Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input_polys_mons,
        const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& result)
{
	SlowMon top;
	for (auto top_var: top_info) {
		top[top_var.index] = top_var.degree;
	}
	auto mon_less = [this](const SlowMon& m1, const SlowMon& m2){return this->MonomialLess(m1, m2);};

	struct MultSlowPolyIterator
	{
		MultSlowPolyIterator(const SlowPolynomial& a_poly)
			:poly(a_poly)
		{}
		SlowMon mul_by;
		SlowPolynomial::const_iterator it;
		SlowPolynomial::const_iterator first_smaller_top_it;//iterator to first element smaller than top (or end)
		const SlowPolynomial& poly;
		SlowMon GetMultMon()
		{
			return Mmul(mul_by, *it);
		}
	};

	struct CombinedFieldFactoryReturningZ2Data
	{
		int max_diferent_numbers_in_coefficients;
		//for more comple creators here will be stored some information about previousely created fields, that i filled on first creation
	};

	struct CombinedFieldFactoryReturningZ2
	{
		ExactFieldAsCombined<ImplementedField> CreateFieldExpectedAsSuitable(CombinedFieldFactoryReturningZ2Data& data, const ImplementedField& field)
		{
			//TODO for non-Z2 case:
			//calculate P_0 based on matrix.size() and max_diferent_numbers_in_coefficients
			//calculate N_0 based on P_0.
			//select field based on P_0 and N_0
			return ExactFieldAsCombined<ImplementedField>(field);
		}
	};

	//prepare:
	//collection of multiplied input polynomials
	//this associates each polynomial in question (corresponding to monomialss of input_polys_mons) with unique number, corresponding to column number
	std::deque<MultSlowPolyIterator> mult_inputs_with_iters;
	//a sorted collection of their monomials in question (only gretear-or-equal than top_info) - each corresponds to matrix row
	std::set<SlowMon, decltype(mon_less)> high_mons(mon_less);
	auto mul_by_it = input_polys_mons.begin();
	auto orig_poly_it = reconstruction_basis.begin();
	for (;orig_poly_it != reconstruction_basis.end(); ++mul_by_it, ++orig_poly_it)
	{
		assert(mul_by_it != input_polys_mons.end());
		for (auto mon: *mul_by_it) {
			auto &mult = emplaced_back(mult_inputs_with_iters, *orig_poly_it);
			SlowMon mul_by;
			for (auto var: mon) {
				mult.mul_by[var.index] = var.degree;
			}
			mult.it = mult.poly.begin();
			for(; mult.it != mult.poly.end(); ++mult.it )
			{
				auto orig_mon = mult.GetMultMon();
				if (mon_less(orig_mon, top))
				{
					break;
				}
				high_mons.insert(orig_mon);
			}
			mult.first_smaller_top_it = mult.it;
		}
	}
	//assert here that top correspond to smallest monomial in high_mons
	assert(!high_mons.empty() && high_mons.find(top) == high_mons.begin());
	ImplementedField field;
	std::vector<std::vector<SparseMatrix::Element<ImplementedField>>> matrix;
	matrix.reserve(high_mons.size());
	//populate matrix rows with (coef from ImplementedField; int column number)
	//move iterator to element smaller than top
	for(auto mult:mult_inputs_with_iters)
	{
		mult.it = mult.first_smaller_top_it;
	}
	for(auto mon:high_mons)
	{
		auto& row =  emplaced_back(matrix);
		for(int col = 0; col < (int)mult_inputs_with_iters.size(); ++col)
		{
			auto mult = mult_inputs_with_iters[col];
			for(;mult.it != mult.poly.begin(); --mult.it)
			{
				auto mult_mon = mult.GetMultMon();
				if (!mon_less(mult_mon, mon))
				{
					if (mon_less(mon, mult_mon))
					{
						break;//too big
					}
					std::remove_reference<decltype(row.back())>::type elem;
					elem.value =  FieldHelpers::One(field);
					elem.column = col;
					row.push_back(elem);
					//eqaul - add to matrix
				}
				//too small find greater
			}
		}
	}
	std::vector<SparseMatrix::Element<ImplementedField>> solution;
	auto factory_data =  Initialized<CombinedFieldFactoryReturningZ2Data>(
			&CombinedFieldFactoryReturningZ2Data::max_diferent_numbers_in_coefficients, std::accumulate(reconstruction_basis.begin(), reconstruction_basis.end(), 0, [](int sum, const SlowPolynomial& poly){return sum + poly.size();})
		);
	//send rows collection to solver that shoud assume that right-side column has 1 in first cell and last rows.size()-1 zeroes
	//solver returns only non-zeros - list of pairs (coef from ImplementedField; int column number)
	SparseMatrix::SolveWithRightSideContainigSingleOne<CombinedFieldFactoryReturningZ2>(field, matrix, solution, factory_data);
	if(solution.empty())
	{
		//if solver fails - return false
		return false;
	}
	auto& new_poly = emplaced_back(*result);
	//add top_info with coef 1.
	new_poly.push_back(top);
	for (auto result_item:solution)
	{
		//calculate sum with coefs given from solver for monomials smaller than top_info, add it to result
		auto& to_mul = mult_inputs_with_iters[result_item.column];
		SlowPolynomial smaller_than_top_part;
		smaller_than_top_part.insert(smaller_than_top_part.begin(), to_mul.first_smaller_top_it, to_mul.poly.end());
		new_poly = PAddZ2(new_poly, smaller_than_top_part, to_mul.mul_by);
	}
	return true;
}

void RingZ2SimpleBase::ConvertResultToFixedMetadataImpl(const unique_deleter_ptr<OutPolysSetForVariyingMetadata>& constructed_result, CrossRingInfo::MonomialListListWithCoef<ImplementedOrder, ImplementedField>& basic_result)
{
	for (auto poly: *constructed_result) {
		basic_result.BeginPolynomialConstruction(poly.size());
		for (auto mon: poly) {
			for (auto var: mon) {
				basic_result.AddVariable(CrossRingInfo::PerVariableData::FromDI(var.second, var.first));
			}
			basic_result.MonomialAdditionDone(FieldHelpers::One(basic_result.Field()));
		}
	}
}

int RingZ2SimpleBase::VarMappingImplReturningOldVarCount(std::vector<int>& new_monomial_vars) const
{

	new_monomial_vars.resize(impl_->keeped_vars_count_ * impl_->new_variables.size());
	int cur_new_var_start = 0;
	for(auto new_var:impl_->new_variables)
	{
		for (auto old_var:new_var)
		{
			new_monomial_vars[cur_new_var_start + old_var.first] = old_var.second;
		}
		cur_new_var_start += impl_->keeped_vars_count_;
	}
	assert(cur_new_var_start == int(new_monomial_vars.size()));
	return impl_->keeped_vars_count_;
}


RingZ2SimpleBase::NewIndices RingZ2SimpleBase::ExtendRingWithMonomialToHelpReconstructImpl(const unique_deleter_ptr<InPolysSetWithOrigMetadata>& reconstruction_basis, Enumerator<CrossRingInfo::PerVariableData> info)
{
	SlowMon& old_mons = emplaced_back(impl_->new_variables);
	for (auto var:info)
	{
		int& new_degree =old_mons[var.index];
		assert(0 == new_degree);
		new_degree = var.degree;
	}
	const int new_var_index = impl_->keeped_vars_count_ + impl_->new_variables.size() - 1;
	SlowMon new_mon;
	new_mon[new_var_index] = 1;
	//new polynomial = old_mons - new_mon
	auto& new_poly = emplaced_back(*reconstruction_basis);
	new_poly.push_back(old_mons);
	new_poly.push_back(new_mon);
	return Initialized<NewIndices>(&NewIndices::new_var_index, new_var_index, &NewIndices::new_poly_index, static_cast<int>(reconstruction_basis->size() - 1));
}

unique_deleter_ptr<RingZ2SimpleBase::OutPolysSetForVariyingMetadata> RingZ2SimpleBase::PrepareEmptyResult()
{
	return MoveToResultType(new OutPolysSetForVariyingMetadata());
}

unique_deleter_ptr<RingZ2SimpleBase::InPolysSetWithOrigMetadata> RingZ2SimpleBase::PrepareForReconstructionImpl(Enumerator<Enumerator<Enumerator<CrossRingInfo::PerVariableData>>> input)
{
	auto result = as_deleter_ptr(new InPolysSetWithOrigMetadata());
	for (auto poly: input) {
		auto& poly_to_populate = emplaced_back(*result);
		for (auto mon: poly) {
			auto& mon_to_populate = emplaced_back(poly_to_populate);
			for (auto var: mon) {
				mon_to_populate[var.index] = var.degree;
			}
		}
	}
	return MoveToResultType(result);
}
