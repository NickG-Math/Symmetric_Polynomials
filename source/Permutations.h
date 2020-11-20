#pragma once
#include "Generators.h"

///@file
///@brief Contains methods for producing permutations, orbits of monomials under symmetric group action, and writing symmetric polynomials in terms of elementary symmetric ones

namespace Symmetric_Polynomials{
///Returns n!
int factorial(int n) {
	if (n == 1)
		return 1;
	else
		return n * factorial(n - 1);
}

///Generates all permutations on a number of letters
class permutations_generator {
	bool completed;
	const char n;
	std::vector<char> generated;
public:
	///Constructs and initializes generator
	permutations_generator(char n) : n(n){
		initialize();
	}
	///Returns 0 if there are no more permutations to be generated
	inline bool keep_going() {
		return (!completed);
	}
	///Initializes generator
	void initialize() {
		generated.resize(n);
		std::iota(generated.begin(), generated.end(), 0);
		completed = 0;
	}
	///Generates new permutation
	void update() {
		if (!std::next_permutation(generated.begin(), generated.end()))
			completed = 1;
	}
	///Returns generated permutation
	std::vector<char> get_generated() const {
		return generated;
	}
};



///Returns vector of all permutations on n letters
std::vector<std::vector<char>> all_permutations(char n) {
	std::vector<std::vector<char>> v;
	v.reserve(factorial(n)); 
	permutations_generator p(n);
	while (p.keep_going()) {
		v.push_back(p.get_generated());
		p.update();
	}
	return v;
}

///Computes the orbit of a monomial under a set of permutations
template<typename scalar_t, typename rel_t>
std::vector<monomial<scalar_t, rel_t>> orbit(const monomial<scalar_t, rel_t>& a, const std::vector<std::vector<char>>& perms) {
	std::vector<monomial<scalar_t, rel_t>> o;
	o.reserve(perms.size());
	for (const auto& i : perms) {
		auto permed = permute(a, i);
		auto it = std::find(o.begin(), o.end(), permed);
		if (it == o.end())
			o.push_back(permed);
	}
	return o;
}

///Returns the greatest monomial in an orbit of permutations
template<typename scalar_t, typename rel_t>
monomial<scalar_t, rel_t> max_in_orbit(const monomial<scalar_t, rel_t>& a, const std::vector<std::vector<char>>& perms) {
	auto max=a;
	for (const auto& i : perms) {
		auto permed = permute(a, i);
		if (max < permed)
			max = permed;
	}
	return max;
}

///Returns a vector of all monomials in given variables and degree organized in orbits under the symmetric group action
template<typename scalar_t, typename rel_t>
std::vector<std::vector<monomial<scalar_t, rel_t>>> all_monomial_orbits(int variables, int degree) {
	auto monos = monomial_basis<scalar_t, rel_t>(variables, degree);
	std::vector<std::vector<monomial<scalar_t,rel_t>>> orbits;
	std::vector<long> already_done;
	orbits.reserve(monos.size());
	already_done.reserve(monos.size());
	orbits.reserve(monos.size());
	std::vector<std::vector<int>> perms=all_permutations(rel_t::true_variables(variables));
	hasher<std::vector<int>> h(std::vector<int>(variables, 0), rel_t::max_powers(variables, degree));
	for (const auto& i : monos) {
		auto hashed = h.hashvector(i.powers);
		auto it = std::find(already_done.begin(), already_done.end(), hashed);
		if (it == already_done.end()) {
			auto o = orbit(i, perms);
			orbits.push_back(o);
			for (const auto& j : o)
				already_done.push_back(h.hashvector(j.powers));
		}
	}
	return orbits;
}

///Class of elementary symmetric polynomials; no relations
template<typename scalar_t>
class elementary_symmetric {
public:
	///Given vector pwr computes e_1^{pwr_1}*...*e_n^{pwr_n}
	polynomial<scalar_t, norelations> compute_product(const std::vector<int>& pwr) {
		polynomial<scalar_t, norelations> product = polynomial<scalar_t, norelations>::constant(1, elementary_symmetric_polynomials.size());
		for (int i = 0; i < pwr.size(); i++)
			if (pwr[i] != 0)
				product = multiply<scalar_t, norelations>(product, power(polynomial<scalar_t, norelations>(orbit(elementary_symmetric_polynomials[i], permutations)), pwr[i], 0));
		return product;
	}

protected:
	std::vector<monomial<scalar_t, norelations>> elementary_symmetric_polynomials; ///<Vector containing the dominant term of the elementary symmetric polynomials e_1,...,e_n (i.e. x_1...x_k)
	std::vector<std::string> elementary_symmetric_names; ///<Vector containing the names "e_i" of the elementary symmetric polynomials e_1,...,e_n

	///Constructs vector of elementary symmetric polynomials in terms of the variables x_1,...,x_n 
	elementary_symmetric(int n) : permutations(all_permutations(n)) {
		elementary_symmetric_polynomials.reserve(n);
		elementary_symmetric_names.reserve(n);
		for (int i = 1; i <= n; i++) {
			std::vector<int> prod_of_first_i(n);
			for (int j = 0; j < i; j++)
				prod_of_first_i[j] = 1;
			elementary_symmetric_polynomials.push_back(monomial<scalar_t, norelations>(1, prod_of_first_i));
			elementary_symmetric_names.push_back("e_" + std::to_string(i));
		}
	}
private:
	const std::vector<std::vector<char>> permutations; ///<All permutations in given number of variables

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///Class for writing a symmetric polynomial as a polynomial on the elementary symmetric polynomials (a "decomposition")
//
///We denote the elementary symmetric polynomials by e_i
//////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename scalar_t>
struct decomposition_elementary_symmetric : public elementary_symmetric<scalar_t> {
	polynomial<scalar_t, norelations> decomposition; ///< The polynomial on the elementary symmetric polynomials.

	///Initialize with number of variables
	decomposition_elementary_symmetric(int n) : elementary_symmetric<scalar_t>(n)  {}

	///Initialize and decompose into polynomial on elementary symmetric polynomials
	decomposition_elementary_symmetric(const polynomial<scalar_t, norelations>& a) : elementary_symmetric<scalar_t>(a.monos[0].powers.size()) { decompose(a);  }


	///Decompose into polynomial on elementary symmetric polynomials
	void decompose(const polynomial<scalar_t, norelations>& a) {
		decomposition.monos.clear();
		decompose_recursive(a);
	}

	///Print a polynomial on elementary symmetric polynomials
	std::string print() const {
		return decomposition.print(this->elementary_symmetric_names);
	}

private:
	void find_power(const monomial<scalar_t, norelations>& term, std::vector<int>& power) {
		for (int i = 0; i < term.powers.size(); i++) {
			if (i == term.powers.size() - 1)
				power.push_back(term.powers[i]);
			else
				power.push_back(term.powers[i] - term.powers[i + 1]);
		}
	}
	void decompose_recursive(const polynomial<scalar_t, norelations>& a) {
		std::vector<int> pwr_here;
		pwr_here.reserve(a.monos[0].powers.size());
		auto max = a.highest_term();
		find_power(max, pwr_here);
		auto product = this->compute_product(pwr_here);
		auto coeff = max.coeff / product.highest_term().coeff;
		decomposition.monos.push_back(monomial<scalar_t, norelations>(coeff,pwr_here));
		auto coeff_product = coeff * product;
		if (a != coeff_product)
			decompose_recursive(a - coeff_product);
	}


};

///Prints decomposition into elementary symmetric polynomials
template<typename scalar_t>
std::ostream& operator<<(std::ostream& os, const decomposition_elementary_symmetric<scalar_t>& a) {
	os << a.print();
	return os;
}
}
