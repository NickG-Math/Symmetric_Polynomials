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
namespace Symmetric_Polynomials {

	///Degree computation given exponent and dimensions (grading)
	template<typename T, typename S>
	int general_compute_degree(const T& exponent, const S& dimensions) {
		int degree = 0;
		for (decltype(exponent.size()) i = 0; i < exponent.size(); i++) {
			degree += exponent[i] * dimensions[i];
		}
		return degree;
	}

	///Returns the sum of all entries of a
	template<typename T, typename S = T>
	S sum(const std::vector<T>& a) {
		S sum=0;
		for (const auto i : a)
			sum += i;
		return sum;
	}

	///Class of rational numbers. 
	//
	///Standard implementation, always keeping numerator and denominator coprime to avoid overflow
	struct Rational {
		long numerator, denominator;
		///Constructs rational n/d in lowest terms
		Rational(long n, long d) {
			auto g = std::gcd(n, d);
			numerator = n / g;
			denominator = d / g;
		}
		///Converts integer to rational
		Rational(long n) : Rational(n, 1) {}
		///Default constructor initializes the 0 rational.
		Rational() : Rational(0) {};
		bool operator ==(const Rational& a) const {
			return (numerator * a.denominator == denominator * a.numerator);
		}
		bool operator !=(const Rational& a) const {
			return !(*this == a);
		}
		Rational operator+(const Rational& a) const {
			return Rational(numerator * a.denominator + denominator * a.numerator, denominator * a.denominator);
		}
		Rational& operator+=(const Rational& a) {
			*this = *this + a;
			return *this;
		}
		Rational operator-(const Rational& a) const {
			return Rational(numerator * a.denominator - denominator * a.numerator, denominator * a.denominator);
		}
		Rational& operator-=(const Rational& a) {
			*this = *this - a;
			return *this;
		}
		Rational operator-() const {
			return Rational(-numerator, denominator);
		}
		Rational operator*(const Rational& a) const {
			return Rational(numerator * a.numerator, denominator * a.denominator);
		}
		Rational& operator*=(const Rational& a) {
			*this = *this * a;
			return *this;
		}
		Rational operator/(const Rational& a) const {
			return Rational(numerator * a.denominator, denominator * a.numerator);
		}
		explicit operator double() const {
			return numerator / denominator;
		}
		bool is_integer() const {
			return !(numerator % denominator);
		}
	};

	///Prints rational number
	std::ostream& operator<<(std::ostream& os, const Rational& a) {
		if (a.is_integer())
			os << a.numerator / a.denominator;
		else
			os << "(" << a.numerator << "/" << a.denominator << ")";
		return os;
	}
}

