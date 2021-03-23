#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#if defined(__GNUC__)
#include "x86intrin.h"
#elif defined(_MSC_VER)
#include "intrin.h"
#endif

///	@file
///	@brief Contains general operations on vectors: hashing, computing degrees

///	@brief The namespace which contains every method and class in the library
namespace symmp
{

	///Hash algorithm CRC32. SSE4.1 and up
	struct crc;

	///Hash algorithm using golden ratio (idea from boost)
	struct boost_hash;

	///	@brief			A generic hashing function that calls other hashing functions.
	///	@tparam	T		The type of element to be hashed (must have a for range loop)
	///	@tparam	hasher	The hashing algorithm. \c boost_hash by default. If SSE4.1 is supported you can also use \c crc
	///	@param	v		The element to be hashed
	///	@return			The hash of the element
	template <typename T, typename hasher = boost_hash>
	size_t generic_hasher(const T &v);

	///	@brief		Degree computation given exponent and dimensions (grading).
	///	@tparam T	The exponent type (eg \c std::vector<int>)
	///	@tparam S	The dimensions type (eg \c std::vector<int>)
	///	@tparam R	The degree type (eg \c uint64_t)
	///	@param exp	The monomial whose degree will be computed, provided via its exponent vector
	///	@param dim	The dimensions of the variables in the monomial
	///	@return		\f$\sum_{i=1}^na_id_i\f$ where exponent=\f$[a_1,...,a_n]\f$ and dimensions=\f$[d_1,...,d_n]\f$
	template <typename R, typename T, typename S>
	R general_compute_degree(const T &exp, const S &dim);

}
#include "impl/General.ipp"