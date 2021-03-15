#pragma once
#include "../General.hpp"
///@file
///@brief Contains general operations on vectors (eg hashing, adding), the class of rational numbers, and methods for applying permutations on vectors


namespace Symmetric_Polynomials {

	//fastest, requires SSE4.1
	struct crc {
		static size_t initialize() { return 0; }
		template<typename T>
		static void combine(size_t& hash, T value) {
			hash = _mm_crc32_u64(hash, value);
		}
	};

	//second fastest
	struct boost_hash {
	public:
		static constexpr uint64_t golden_ratio = 0x9e3779b97f4a7c15;
		static size_t initialize() { return 0; }
		template<typename T>
		static void combine(size_t& hash, T value) {
			hash ^= value + golden_ratio + (hash << 6) + (hash >> 2);
		}
	};

	template<typename T, typename hasher>
	size_t generic_hasher(const T& v) {
		auto hash = hasher::initialize();
		for (const auto i : v)
			hasher::combine(hash, i);
		return hash;
	}

	template<typename R, typename T, typename S>
	R general_compute_degree(const T& exponent, const S& _dimensions) {
		R degree = 0;
		for (size_t i = 0; i < exponent.size(); i++)
			degree += exponent[i] * _dimensions[i];
		return degree;
	}

}

