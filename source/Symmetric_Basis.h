#pragma once
#include "Polynomials.h"
#include "Generators.h"

///@file
///@brief Contains methods for writing symmetric polynomials (with no relations) in terms of the elementary symmetric ones. Also contains a factory class that generalizes the interface to any subring of a polynomial ring

namespace Symmetric_Polynomials {

	///The standard variables \f$x_i\f$ in a polynomial, with \f$|x_i|=1\f$ and no relations.
	template<typename T = int64_t, typename _deg = int64_t>
	struct Standard_Variables : public std::vector<T> {
		using std::vector<T>::vector;
		typedef _deg deg_t; ///<Degree typedef
		///Degree of exponent
		deg_t degree() const {
			deg_t sum = 0;
			for (const auto i : *this)
				sum += i;
			return sum;
		}
		///Name of variable
		static std::string name(int i, int number_of_variables) {
			return "x_" + std::to_string(i + 1);
		}
		///Addition of exponents
		Standard_Variables operator+ (const Standard_Variables& other) const {
			Standard_Variables v(*this);
			for (size_t i = 0; i < this->size(); i++)
				v[i] += other[i];
			return v;
		}
		///Standard hash function
		size_t operator ()() const {
			return generic_hasher(*this);
		}
	};

	///Variables \f$e_1,...,e_n\f$ denoting the elementary symmetric polynomials \f$e_i=\sigma_i\f$ of degrees \f$|e_i|=i\f$
	template<typename T = int64_t, typename _deg = int64_t>
	struct Elementary_Symmetric_Variables : public Standard_Variables<T, _deg> {
		using Standard_Variables<T, _deg>::Standard_Variables;
		typedef _deg deg_t; ///<Degree typedef
		///The degree of a monomial on the \f$e_i\f$, \f$|e_i|=i\f$
		deg_t degree() const {
			deg_t deg = 0;
			for (size_t i = 0; i < this->size(); i++)
				deg += (*this)[i] * (i + 1);
			return deg;
		}
		///The names of the \f$e_i\f$
		static std::string name(int i, int number_of_variables) {
			return "e_" + std::to_string(i + 1);
		}
		///Standard hash function
		size_t operator ()() const {
			return generic_hasher(*this);
		}
	};


	///////////////////////////////////////////////////////////////////////////////////
	///Factory class that provides the general interface of a generating basis for a subring of a polynomial ring. Example implementations are Symmetric_Basis and Half_Idempotent_Basis.
	//
	///Use by inheriting and setting T be the inheriting class (CRTP). Then construct _generators (an optionally generator_names and generator_dimensions) in the inheriting class.
	///
	///The inheriting class T must implement a method find_exponent with singature: Polynomial<container_2_t>::exp_t find_exponent(const Polynomial<container_1_t>::exp_t&);
	///
	///container_1_t is the container type of the original variables (of the polynomial ring)
	///
	///container_2_t is the container type of the new variables that are the basis elements (of the subring)
	///////////////////////////////////////////////////////////////////////////////////
	template<typename T, typename container_1_t, typename container_2_t>
	class Polynomial_Basis {
	public:

		///Transform a polynomial on the original variables to one on the generating basis
		Polynomial<container_2_t> operator()(Polynomial<container_1_t> a) const {

			//set dimensions and names if nonempty
			const typename Polynomial<container_2_t>::deg_t* gen_dims = generator_dimensions.empty() ? nullptr : generator_dimensions.data();
			const std::string* gen_names = generator_names.empty() ? nullptr : generator_names.data();
			Polynomial<container_2_t> decomposition(number_of_variables, gen_dims, gen_names);

			while (true) {
				auto max = a.highest_term();
				auto exponent = static_cast<const T*>(this)->find_exponent(max.exponent());
				auto product = compute_product(exponent);
				auto coeff = max.coeff() / product.highest_term().coeff();
				decomposition.insert(exponent, coeff);
				product *= coeff;
				if (a == product)
					return decomposition;
				else
					a -= product;
			}
		}

		///Transform a polynomial on the generating basis into a polynomial on the original variables
		Polynomial<container_1_t> operator()(const Polynomial<container_2_t>& a) const {
			Polynomial<container_1_t> p(number_of_variables);
			for (auto it = a.begin(); it != a.end(); ++it) {
				auto prod = compute_product(it.exponent());
				prod *= it.coeff();
				p += prod;
			}
			return p;
		}

		///Constructor given number of variables
		Polynomial_Basis(int number_of_variables) : number_of_variables(number_of_variables) {}

		///Returns const& of vector containing the generating basis
		const auto& generators() const {
			return _generators;
		}

		///Returns const& of vector containing the dimensions of the generating basis (can be empty!)
		const auto& dimensions() const {
			return generator_dimensions;
		}

		///Returns const& of vector containing the names of the generating basis (can be empty!)
		const auto& names() const {
			return generator_names;
		}

		///The number of (the original) variables of the polynomial ring
		const int number_of_variables;

	protected:

		///The generators of the polynomial basis, constructed in the inheriting class
		std::vector<Polynomial<container_1_t>> _generators;

		///The dimensions of the generators, optionally constructed in the inheriting class
		std::vector<typename Polynomial<container_2_t>::deg_t> generator_dimensions;

		///The names of the generators, optionally constructed in the inheriting class
		std::vector<std::string> generator_names;

	private:

		Polynomial<container_1_t> compute_product(const typename Polynomial<container_2_t>::exp_t& exponent) const {
			Polynomial<container_1_t> product(number_of_variables, 1);
			for (size_t i = 0; i < _generators.size(); i++)
				if (exponent[i] != 0)
					product *= (_generators[i] ^ exponent[i]);
			return product;
		}
	};

	///Class for symmetric polynomials with no relations, allowing transformation from \f$x_i\f$ variables to \f$e_i\f$ variables and vice-versa.
	template<typename, typename>
	class Symmetric_Basis;

	template<typename x_container_t, typename e_container_t>
	class Symmetric_Basis : public Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t> {

		using Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t>::number_of_variables;
		using Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t>::_generators;

	public:
		///Constructor given number of variables
		Symmetric_Basis(int n) : Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t>(n) {
			_generators.reserve(number_of_variables);
			for (int i = 1; i <= number_of_variables; i++)
				_generators.push_back(get_elementary_symmetric(i));	
		}

	private:

		typedef typename Polynomial<x_container_t>::exp_t x_t;
		typedef typename Polynomial<e_container_t>::exp_t e_t;

		Polynomial<x_container_t> get_elementary_symmetric(int i) const {
			Polynomial<x_container_t> poly(number_of_variables);
			x_t mono(number_of_variables);
			auto c = Combination_Generator(number_of_variables, i);
			for (const auto& comb : c) {
				for (const auto j : comb)
					mono[j] = 1;
				poly.insert(mono, 1);
				for (const auto j : comb)
					mono[j] = 0;
			}
			return poly;
		}

		e_t find_exponent(const x_t& term) const {
			e_t exponent;
			exponent.reserve(number_of_variables);
			for (int i = 0; i < term.size(); i++) {
				if (i == term.size() - 1)
					exponent.push_back(term[i]);
				else
					exponent.push_back(term[i] - term[i + 1]);
			}
			return exponent;
		}

		///Befriending parent for CRTP.
		friend class Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t>;
	};


	///Aliases Symmetric_Basis for certain default template parameters
	template<typename scl_t, typename exp_value_t = int64_t, typename deg_t = int64_t, bool ordered = 1>
	using Symmetric_Basis_Default = Symmetric_Basis<default_container<scl_t, Standard_Variables<exp_value_t, deg_t>, ordered>, default_container<scl_t, Elementary_Symmetric_Variables<exp_value_t, deg_t>, ordered>>;

}
