#pragma once
#include "General.h"
#include <functional>
///@file
///@brief Contains the classes of monomials and polynomials in multiple variables with relations. 

namespace Symmetric_Polynomials{
int general_compute_degree(const std::vector<int>& powers, const std::vector<int>& dimensions) {
	int degree = 0;
	for (int i = 0; i < powers.size(); i++) {
		degree += powers[i] * dimensions[i];
	}
	return degree;
}

template<typename scalar_t, typename rel_t>
class polynomial;

///////////////////////////////////////////////////////////////////////////////////
///Class of monomials in multiple variables with relations.
//
///A monomial \f$a=cx_1^{p_1}x_2^{p_2}....x_n^{p_n}\f$ is defined by its coefficient \f$c\f$ and powers vector \f$[p_1,...,p_n]\f$.
///The template scalar_t determines the type of the coefficient \f$c\f$.
///The template rel_t determines the relations that the variables \f$x_i\f$ satisfy.
///////////////////////////////////////////////////////////////////////////////////////
template<typename scalar_t, typename rel_t>
class monomial {
private:
	void getdegree() {
		if (coeff != 0)
			degree = rel_t::compute_degree(powers);
		else
			degree = 0;
	}
	int degree; ///<Degree of monomial

public:

	scalar_t coeff; ///<Coefficient of monomial
	std::vector<int> powers; ///<Powers vector

	///Constructs monomial given coefficient and powers
	monomial(const scalar_t& coeff, const std::vector<int>& powers) :
		coeff(coeff), powers(powers) {
		rel_t::apply_monomial(this);
		getdegree();
	};

	///Default constructor
	monomial() {};


	///Constructs constant monomial with given coefficient and number of variables
	static monomial<scalar_t, rel_t> constant(int coeff, int variables) {
		return monomial<scalar_t, rel_t>(coeff, std::vector<int>(variables)); //monomial 1
	}

	///Returns product of monomials; depending on the relations, this could be a polynomial after applying them, so the return type is conditional on that
	template<typename s_scalar_t, typename s_rel_t>
	inline std::conditional_t<s_rel_t::product_of_monomials_is_monomial, monomial<s_scalar_t, s_rel_t>, polynomial<s_scalar_t, s_rel_t>> operator *(const monomial<s_scalar_t, s_rel_t>& b) const {
		return monomial<s_scalar_t, s_rel_t>(this->coeff * b.coeff, this->powers + b.powers);
	}

	///Returns -b for monomial b
	monomial<scalar_t, rel_t> operator -() const {
		auto b = *this;
		b.coeff = -b.coeff;
		return b;
	}

	///Division of monomials. No checks are made to ensure that the resulting powers are nonnegative
	monomial<scalar_t, rel_t> operator /(const monomial<scalar_t, rel_t>& b) const {
		return monomial<scalar_t, rel_t>(this->coeff / b.coeff, this->powers - b.powers);
	}

	///Standard equality of monomials
	bool operator==(const monomial<scalar_t, rel_t>& b) const {
		return (this->coeff == b.coeff && this->powers == b.powers);
	}

	///a<b iff degree(a)<degree(b) or they have the same degree and for the powers vectors, powers(a)<powers(b)
	bool operator<(const monomial<scalar_t, rel_t>& b) const {
		return (this->degree < b.degree) || (this->degree == b.degree && this->powers < b.powers);
	}

	///Print monomial using given variable names
	std::string print(const std::function<std::string(int)>& variable_names) const {
		bool started = 0;
		std::stringstream ss;
		if (coeff != 1)
			ss << coeff << "*";
		for (int i = 0; i < powers.size(); i++) {
			if (powers[i]>1)
				ss << variable_names(i) << "^" << powers[i] << "*";
			else if (powers[i]==1)
				ss << variable_names(i) << "*";
		}
		auto str = ss.str();
		if (str.back() == '*')
			str.pop_back();
		return str;
	}

	///Print monomial using given variable names
	std::string print(const std::vector<std::string>& variable_names) const {
		std::function<std::string(int)> f = [=](int i) {return variable_names[i];};
		return print(f);
	}

	///Print monomial using default variable names x_1,...,x_n
	std::string print() const {
		std::function<std::string(int)> f = [](int i) {return "x_" + std::to_string(i + 1);};
		return print(f);
	}

	template<typename s_scalar_t, typename s_rel_t>
	friend class polynomial;

};

///Product of scalar and monomial
template<typename scalar_t, typename rel_t>
inline monomial<scalar_t, rel_t> operator *(const scalar_t& a, const monomial<scalar_t, rel_t>& b) {
	return monomial<scalar_t, rel_t>(a * b.coeff, b.powers);
}

///Comparator for monomials based on <
template<typename scalar_t, typename rel_t>
struct monomial_compare
{
	inline bool operator() (const monomial<scalar_t, rel_t>& a, const monomial<scalar_t, rel_t>& b) {
		return (a < b);
	}
};

///Prints monomial
template<typename scalar_t, typename rel_t>
std::ostream& operator<<(std::ostream& os, const monomial<scalar_t, rel_t>& a) {
	os << a.print();
	return os;
}

///Permutes the variables in a monomial according permutation; returns permuted monomial
template<typename scalar_t, typename rel_t>
monomial<scalar_t, rel_t> permute(const monomial<scalar_t, rel_t>& a, const std::vector<char>& perm) {
	return monomial<scalar_t, rel_t>(a.coeff, rel_t::permute(a.powers, perm));
}

/////////////////////////////////////////////////////////////////////////////
///Class of polynomials in multiple variables with relations.
//
/// A polynomial is just a vector of monomials, ordered in increasing degree
//
///The template scalar_t determines the type of the coefficients.
///The template rel_t determines the relations that the variables x_i satisfy.
//////////////////////////////////////////////////////////////////////////////////////
template<typename scalar_t, typename rel_t>
class polynomial {
	void reorder() {
		std::sort(monos.begin(), monos.end(), monomial_compare<scalar_t, rel_t>());
	}

public:
	///The vector of monomials, ordered in increasing degree
	std::vector<monomial<scalar_t, rel_t>> monos; 

	///Incorporates all monomials with same powers vector into a single monomial by adding coefficients, condensing the polynomial
	void condense() {
		int i = 0;
		while (i < monos.size()) {
			if (monos[i].coeff == 0) //erase 0 monomial and dont change i
				monos.erase(monos.begin() + i);
			else if (i != monos.size() - 1 && monos[i].degree == monos[i + 1].degree && monos[i].powers == monos[i + 1].powers) {
				if (monos[i].coeff + monos[i + 1].coeff == 0) { //they cancel out so erase both and dont change i
					monos.erase(monos.begin() + i, monos.begin() + i + 2);
				}
				else {
					monos[i].coeff += monos[i + 1].coeff; //incorporate to the first and erase the second and dont change i
					monos.erase(monos.begin() + i + 1);
				}
			}
			else { //nothing erased, change i
				i++;
			}
		}
	}

	///Returns highest term monomial in a polynomial
	monomial<scalar_t, rel_t> highest_term() const {
		return monos.back();
	}

	///Constructs a polynomial given a vector of monomials, applying relations to variables and reordering as needed. Condensation is performed by default but can be disabled at the user's risk
	polynomial(const std::vector<monomial<scalar_t, rel_t>>& monos, bool dontcondense = 0) : monos(monos) {
		rel_t::apply_polynomial(this);
		reorder();
		if (!dontcondense)
			condense();
	}

	///Constructs a polynomial from a single monomial
	polynomial(const monomial<scalar_t, rel_t>& mono) : polynomial<scalar_t, rel_t>(std::vector<monomial<scalar_t, rel_t>>({ mono })) {}

	///Constructs constant polynomial given coefficient and number of variables
	static polynomial<scalar_t, rel_t> constant(int coeff, int variables) {
		return polynomial<scalar_t, rel_t>(monomial<scalar_t, rel_t>::constant(coeff, variables)); //polynomial 1
	}


	///Default Constructor
	polynomial() {}


	///Adds two polynomials
	polynomial<scalar_t, rel_t> operator+(const polynomial<scalar_t, rel_t>& b) const {
		auto sum = join(monos, b.monos);
		return polynomial<scalar_t, rel_t>(sum);
	}

	///Returns -b for polynomial b
	polynomial<scalar_t, rel_t> operator-() const {
		auto b = *this;
		for (auto& i : b.monos)
			i.coeff = -i.coeff;
		return b;
	}

	///Subtracts polynomial from polynomial
	polynomial<scalar_t, rel_t> operator-(const polynomial<scalar_t, rel_t>& b) const {
		return *this + (-b);
	}

	///Standard equality of polynomials
	bool operator==(const polynomial<scalar_t, rel_t>& b) const {
		return (monos == b.monos);
	}

	///Standard inequality of polynomials
	bool operator!=(const polynomial<scalar_t, rel_t>& b) const {
		return !(*this == b);
	}


	///Print polynomial using given variable names
	std::string print(const std::function<std::string(int)>& variable_names) const {
		std::stringstream ss;
		for (int i = 0; i < monos.size(); i++)
			if (i != 0)
				ss << " + " << monos[i].print(variable_names);
			else
				ss << monos[i].print(variable_names);
		return ss.str();
	}

	///Print monomial using given variable names
	std::string print(const std::vector<std::string>& variable_names) const {
		std::function<std::string(int)> f = [=](int i) {return variable_names[i];};
		return print(f);
	}

	///Print polynomial using default variable names x_1,...,x_n
	std::string print() const {
		std::stringstream ss;
		for (int i = 0; i < monos.size(); i++)
			if (i != 0)
				ss << " + " << monos[i].print();
			else
				ss << monos[i].print();
		return ss.str();
	}
};

///Multiplies two polynomials with the ability to disable condensation (not recommended)
template<typename scalar_t, typename rel_t>
polynomial<scalar_t, rel_t> multiply(const polynomial<scalar_t, rel_t>& a, const polynomial<scalar_t, rel_t>& b, bool dontcondense = 0) {
	if (a.monos.empty()) //0 polynomial
		return a;
	if (b.monos.empty()) //0 polynomial
		return b;
	std::vector<monomial<scalar_t, rel_t>> prodmonos;
	prodmonos.reserve(a.monos.size() * b.monos.size());
	for (const auto& i : a.monos)
		for (const auto& j : b.monos)
			if (i.coeff != 0 && j.coeff != 0) {
				if constexpr (rel_t::product_of_monomials_is_monomial)
					prodmonos.push_back(i * j);
				else {
					polynomial<scalar_t, rel_t> ij = i * j;
					for (const auto& k : ij.monos)
						prodmonos.push_back(k);
				}
			}
	return polynomial<scalar_t, rel_t>(prodmonos, dontcondense);
}

///Multiplies scalar and polynomial
template<typename scalar_t, typename rel_t>
polynomial<scalar_t, rel_t> operator*(const scalar_t& a, const polynomial<scalar_t, rel_t>& b) {
	if (a == 0) //0 polynomial
		return polynomial<scalar_t, rel_t>::constant(0, b.monos[0].powers.size());
	if (a == 1)
		return b;
	auto ab = b;
	for (auto& i : ab.monos)
		i.coeff *= a;
	return ab;
}

template<typename scalar_t, typename rel_t>
inline polynomial<scalar_t, rel_t> operator*(const polynomial<scalar_t, rel_t>& a, const polynomial<scalar_t, rel_t>& b) {
	return multiply(a, b);
}

///Returns the product of monomial and polynomial
template<typename scalar_t, typename rel_t>
inline polynomial<scalar_t, rel_t> operator*(const monomial<scalar_t, rel_t>& a, const polynomial<scalar_t, rel_t>& b) {
	return multiply(polynomial<scalar_t, rel_t>(a), b);
}

///Returns the product of polynomial and monomial
template<typename scalar_t, typename rel_t>
inline polynomial<scalar_t, rel_t> operator*(const polynomial<scalar_t, rel_t>& a, const monomial<scalar_t, rel_t>& b) {
	return multiply(a, polynomial<scalar_t, rel_t>(b));
}

///Computes a^p for polynomial a and nonnegative integer p
template<typename scalar_t, typename rel_t>
polynomial<scalar_t, rel_t> power(const polynomial<scalar_t, rel_t>& a, int p, bool dontcondense = 0) {
	if (p == 0)
		return polynomial<scalar_t, rel_t>::constant(1, a.monos[0].powers.size()); //polynomial 1
	auto pwr = a;
	for (int i = 2; i <= p;i++)
		pwr = multiply(pwr, a, dontcondense);
	return pwr;
}

///Prints polynomial
template<typename scalar_t, typename rel_t>
std::ostream& operator<<(std::ostream& os, const polynomial<scalar_t, rel_t>& a) {
	os << a.print();
	return os;
}
}
