#pragma once
#include "Polynomials.h"

///@file
///@brief Contains two classes of relations on polynomial variables: no relations, and when half the variables are idempotents. Possible to extend by adding new types of relations

namespace Symmetric_Polynomials{
///Base class for all relations. Inherit from it and specialize static methods to get a functioning relations class
struct relations_base {

	///By default 0, but set to 1 in a specialization if you can guarantee that the product of monomials after relations are applied is a monomial (as opposed to polynomial)
	inline static const bool product_of_monomials_is_monomial=0;

	///Compute degree based on the powers vector and relations
	static int compute_degree(const std::vector<int>&) { return -1; }

	///Apply relations on monomial
	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr) {};

	///Apply relations on polynomial
	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	///Return max powers vector corresponding to polynomial of given degree and number of variables
	static std::vector<int> max_powers(int, int);

	///Apply permutation on powers vector
	static std::vector<int> permute(const std::vector<int>&, const std::vector<char>&);

	///Actual number of variables given length of powers vector
	static int true_variables(int n);
};


///Class when no relations are needed
struct norelations : relations_base {

	inline static const bool product_of_monomials_is_monomial = 1;

	static int compute_degree(const std::vector<int>& powers) {
		return sum(powers);
	}

	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr) {}

	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	static std::vector<int> max_powers(int variables, int degree) {
		return std::vector<int>(variables, degree);
	}

	static std::vector<int> permute(const std::vector<int>& v, const std::vector<char>& perm) {
		return apply_permutation(v, perm);
	}

	static int true_variables(int n) { return n; }
};

///Class when we have variables \f$x_1,...,x_n,y_1,...,y_n\f$ and the relations are: \f$y_i^2=y_i\f$
struct halfidempotent : relations_base {
	inline static const bool product_of_monomials_is_monomial = 1;

	static int compute_degree(const std::vector<int>& powers) {
		return sum(powers,powers.size()/2);
	}

	template<typename monomial_ptr>
	static void apply_monomial(monomial_ptr m) {
		for (int i = m->powers.size() / 2; i < m->powers.size(); i++)
			m->powers[i] = std::min(m->powers[i], 1);
	}

	template<typename polynomial_ptr>
	static void apply_polynomial(polynomial_ptr) {};

	static std::vector<int> max_powers(int variables, int degree) {
		std::vector<int> max(variables, degree);
		for (int i = max.size() / 2; i < max.size(); i++)
			max[i] = 1;
		return max;
	}
	static std::vector<int> permute(const std::vector<int>& v, const std::vector<char>& perm) {
		return apply_permutation_pieces(v, 2, perm);
	}
	static int true_variables(int n) { return n / 2; }
};
}
