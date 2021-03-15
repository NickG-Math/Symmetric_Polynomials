#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#if defined (__GNUC__)
#include "x86intrin.h"
#elif defined(_MSC_VER)
#include "intrin.h"
#endif

///@file
///@brief Contains general operations on vectors: hashing, adding

///The namespace which contains every method and class in the library
namespace Symmetric_Polynomials {

	struct crc; ///<Hash algorithm CRC32
	struct boost_hash; ///<Hash algorithm using golden ratio (idea from boost)

	///A generic hashing function that calls other hashing functions.
	//
	///@tparam T The type of element to be hashed (must have a for range loop)
	///@tparam hasher The hasing algorithm. boost_hash by default
	///@param v The element to be hashed
	///@return The hash of the element
	template<typename T, typename hasher= boost_hash>
	size_t generic_hasher(const T& v);


	///Degree computation given exponent and dimensions (grading). 
	//
	///@tparam T The exponent type (eg std::vector<int>)
	///@tparam S The dimensions type (eg std::vector<int>)
	///@tparam R The degree type (eg uint64_t)
	///@param exponent The monomial whose degree will be computed, provided via its exponent vector
	///@param dimensions The dimensions of the variables in the monomial
	///@return \f$\sum_{i=1}^na_id_i\f$ where exponent=\f$[a_1,...,a_n]\f$ and dimensions=\f$[d_1,...,d_n]\f$
	template<typename R, typename T, typename S>
	R general_compute_degree(const T& exponent, const S& dimensions);

}
#include "impl/General.ipp"
