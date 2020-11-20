#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#include <algorithm>
///@file
///@brief Contains general operations on vectors (eg hashing, adding), the class of rational numbers, and methods for applying permutations on vectors

///The namespace which contains every method and class in the library
namespace Symmetric_Polynomials{

///Hashes all vectors between two vectors min and max into a long integer
template<typename T>
struct hasher {
	T min, max_minus_min_plus_1;
	hasher(const T& min, const T& max) : min(min) {
		max_minus_min_plus_1.reserve(min.size());
		for (int i = 0; i < min.size(); i++)
			max_minus_min_plus_1.push_back(max[i] - min[i] + 1);
	};
	long hashvector(const T& deg) {
		long hash = deg[0] - min[0];
		long prod = max_minus_min_plus_1[0];
		for (long i = 1; i < deg.size(); i++) {
			hash += (deg[i] - min[i]) * prod;
			prod *= max_minus_min_plus_1[i];
		}
		return hash;
	}
};


///Applies permutation perm on the target from starting point begin
template<typename Iterator>
void apply_permutation(std::vector<int>& target, Iterator begin, const std::vector<char>& perm) {
	for (const auto i : perm)
		target.push_back(*(begin + i));
}


///Breaks vector v into given number of pieces and applies given permutation on each one; then joins the pieces and returns the vector
std::vector<int> apply_permutation_pieces(const std::vector<int>& v, int pieces, const std::vector<char>& perm) {
	std::vector<int> target;
	target.reserve(v.size());
	for (int i = 0; i < pieces; i++)
		apply_permutation(target, v.begin() + i * perm.size(), perm);
	return target;
}

///Applies permutation perm on the vector v
std::vector<int> apply_permutation(const std::vector<int>& v, const std::vector<char>& perm) {
	return apply_permutation_pieces(v, 1, perm);
}


///For two vectors a,b,  a<=b iff a[i]<=b[i] for all i
template<typename T, typename S>
bool operator <=(const std::vector<T>& a, const std::vector<S>& b) {
	for (int i = 0; i < a.size(); i++)
		if (a[i] > b[i])
			return 0;
	return 1;
}

///Returns a[0]+...+a[end]
template<typename T, typename S = T>
S sum(const std::vector<T>& a, int end) {
	S sum = a[0];
	for (int i = 1; i < end; i++)
		sum = sum + S(a[i]);
	return sum;
}

///Returns the sum of all entries of a
template<typename T, typename S = T>
S sum(const std::vector<T>& a) {
	return sum<T, S>(a, a.size());
}


///Adds two vectors coordinate wise
std::vector<int> operator +(const std::vector<int>& a, const std::vector<int>& b) {
	std::vector<int> sum;
	sum.reserve(a.size());
	for (int i = 0; i < a.size(); i++)
		sum.push_back(a[i] + b[i]);
	return sum;
}

///Subtract a vector from another coordinate wise
std::vector<int> operator -(const std::vector<int>& a, const std::vector<int>& b) {
	std::vector<int> diff;
	diff.reserve(a.size());
	for (int i = 0; i < a.size(); i++)
		diff.push_back(a[i] - b[i]);
	return diff;
}
 
///Multiplies scalar and vector coordinate wise
std::vector<int> operator *(int a, const std::vector<int>& b) {
	std::vector<int> prod;
	prod.reserve(b.size());
	for (const auto& i : b)
		prod.push_back(a * i);
	return prod;
}

///Joins two vectors a,b into a new vector with entries a[0],...,a[a.size()-1],b[0],...,b[b.size()-1]
template<typename T>
std::vector<T> join(const std::vector<T>& a, const std::vector<T>& b) {
	std::vector<T> joined;
	joined.reserve(a.size() + b.size());
	joined.insert(joined.end(), a.begin(), a.end());
	joined.insert(joined.end(), b.begin(), b.end());
	return joined;
}


///Class of rational numbers. Standard implementation
struct rational {
	int numerator, denominator;
	///Constructs rational n/d in lowest terms
	rational(int n, int d) {
		auto g = std::gcd(n, d);
		numerator = n / g;
		denominator = d / g;
	}
	///Converts integer to rational
	rational(int n) : rational(n, 1) {}
	///Default constructor initializes the 0 rational.
	rational() : rational(0) {};
	bool operator ==(const rational& a) const {
		return (numerator * a.denominator == denominator * a.numerator);
	}
	bool operator !=(const rational& a) const {
		return !(*this == a);
	}
	rational operator+(const rational& a) const {
		return rational(numerator * a.denominator + denominator * a.numerator, denominator * a.denominator);
	}
	rational& operator+=(const rational& a) {
		*this = *this + a;
		return *this;
	}
	rational operator-(const rational& a) const {
		return rational(numerator * a.denominator - denominator * a.numerator, denominator * a.denominator);
	}
	rational operator-() const {
		return rational(-numerator, denominator);
	}
	rational operator*(const rational& a) const {
		return rational(numerator * a.numerator, denominator * a.denominator);
	}
	rational& operator*=(const rational& a) {
		*this = *this * a;
		return *this;
	}
	rational operator/(const rational& a) const {
		return rational(numerator * a.denominator, denominator * a.numerator);
	}
	explicit operator double() const {
		return numerator / denominator;
	}
	bool is_integer() const {
		return !(numerator % denominator);
	}
};

///Prints rational number
std::ostream& operator<<(std::ostream& os, const rational& a) {
	if (a.is_integer())
		os << a.numerator / a.denominator;
	else
		os << "(" << a.numerator << "/" << a.denominator << ")";
	return os;
}
}

