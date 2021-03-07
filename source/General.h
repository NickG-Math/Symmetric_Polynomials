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
///@brief Contains general operations on vectors (eg hashing, adding), the class of rational numbers, and methods for applying permutations on vectors


namespace { //some hashing functions

	//fastest, requires SSE4.1
	struct crc {
		static size_t initialize() { return 0; }
		template<typename T>
		static void combine(size_t& hash, T value) {
			hash = _mm_crc32_u64(hash, value);
		}
	}; 

	//second fastest
	class boost_hash {
		static constexpr uint64_t golden_ratio = 0x9e3779b97f4a7c15;
	public:
		static size_t initialize() { return 0; }
		template<typename T>
		static void combine(size_t& hash, T value) {
			hash ^= value + golden_ratio + (hash << 6) + (hash >> 2);
		}
	};

}


///The namespace which contains every method and class in the library
namespace Symmetric_Polynomials {

	///A generic hashing function that calls other hashing functions (second template parameter)
	template<typename T, typename hasher= boost_hash>
	size_t generic_hasher(const T& v) {
		auto hash=hasher::initialize();
		for (const auto i:v)
			hasher::combine(hash, i);
		return hash;
	}


	///Degree computation given exponent and dimensions (grading)
	template<typename R, typename T, typename S>
	R general_compute_degree(const T& exponent, const S& _dimensions) {
		R degree = 0;
		for (size_t i = 0; i < exponent.size(); i++)
			degree += exponent[i] * _dimensions[i];
		return degree;
	}

}

