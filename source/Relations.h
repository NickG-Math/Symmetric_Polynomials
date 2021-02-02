#pragma once
#include "Polynomials.h"

///@file
///@brief Contains two classes of relations on polynomial variables: no relations, and when half the variables are idempotents. Possible to extend by adding new types of relations

namespace Symmetric_Polynomials{
///Base class for all relations. Inherit from it and specialize static methods to get a functioning relations class
struct relations_base {

	///By default 0, but set to 1 in a specialization if you can guarantee that the product of monomials after relations are applied is a monomial (as opposed to polynomial)
	inline static const bool product_of_monomials_is_monomial=0;

	///Compute degree based on the exponent vector and relations
	static int compute_degree(const std::vector<char>&) { return -1; }

	///Apply relations on monomial
	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr) {};

	///Apply relations on polynomial
	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	///Return max exponent vector corresponding to polynomial of given degree and number of variables
	static std::vector<char> max_exponent(int, int);

	///Apply permutation on exponent vector
	static std::vector<char> permute(const std::vector<char>&, const std::vector<char>&);

	///Actual number of variables given length of exponent vector
	static int true_variables(int n);
};


///Class when no relations are needed
struct norelations : relations_base {

	///Product of monomials with no relations is a monomial
	inline static const bool product_of_monomials_is_monomial = 1;
	
	///The degree of a monomial with given exponent
	static int compute_degree(const std::vector<char>& exponent) {
		return sum(exponent);
	}

	///Apply relations on monomial (i.e. do nothing)
	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr) {}

	///Apply relations on polynomial  (i.e. do nothing)
	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	///Return max exponent vector corresponding to polynomial of given degree and number of variables
	static std::vector<char> max_exponent(int variables, int degree) {
		return std::vector<char>(variables, degree);
	}

	///Apply permutation on exponent vector
	static std::vector<char> permute(const std::vector<char>& v, const std::vector<char>& perm) {
		return apply_permutation(v, perm);
	}

	///Actual number of variables given length of exponent vector (i.e. x_1,...,x_n and n is that number)
	static int true_variables(int n) { return n; }
};

///Class when we have variables \f$x_1,...,x_n,y_1,...,y_n\f$ and the relations are: \f$y_i^2=y_i\f$
struct halfidempotent : relations_base {
	///Product of monomials with half idempotent relations is a monomial
	inline static const bool product_of_monomials_is_monomial = 1;

	///The degree of a monomial with given exponent
	static int compute_degree(const std::vector<char>& exponent) {
		return sum(exponent,exponent.size()/2);
	}

	///Apply relations on monomial
	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr m) {
		for (int i = m->exponent.size() / 2; i < m->exponent.size(); i++)
			m->exponent[i] = std::min(m->exponent[i], (char)1);
	}

	///Apply relations on polynomial
	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	///Return max exponent vector corresponding to polynomial of given degree and number of variables
	static std::vector<char> max_exponent(int variables, int degree) {
		std::vector<char> max(variables, degree);
		for (int i = max.size() / 2; i < max.size(); i++)
			max[i] = 1;
		return max;
	}
	///Apply permutation on exponent vector, separately on the variables x_i and the variables y_i
	static std::vector<char> permute(const std::vector<char>& v, const std::vector<char>& perm) {
		return apply_permutation_pieces(v, 2, perm);
	}
	///Actual number of variables given length of exponent vector (i.e. x_1,...,x_n,y_1,...,y_n and n is that number)
	static int true_variables(int n) { return n / 2; }
};
}
