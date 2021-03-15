#pragma once
#include "../Symmetric_Basis.hpp"

///@file
///@brief Contains methods for writing symmetric polynomials (with no relations) in terms of the elementary symmetric ones. Also contains a factory class that generalizes the interface to any subring of a polynomial ring

namespace Symmetric_Polynomials {


	template<typename T, typename _deg>
	_deg Standard_Variables<T, _deg>::degree() const {
		deg_t sum = 0;
		for (const auto i : *this)
			sum += i;
		return sum;
	}

	template<typename T, typename _deg>
	std::string Standard_Variables<T, _deg>::name(int i, int number_of_variables) {
		return "x_" + std::to_string(i + 1);
	}

	template<typename T, typename _deg>
	Standard_Variables<T, _deg> Standard_Variables<T, _deg>::operator+ (const Standard_Variables& other) const {
		Standard_Variables v(*this);
		for (size_t i = 0; i < this->size(); i++)
			v[i] += other[i];
		return v;
	}

	template<typename T, typename _deg>
	size_t Standard_Variables<T, _deg>::operator ()() const {
		return generic_hasher(*this);
	}

	template<typename T, typename _deg>
	_deg Elementary_Symmetric_Variables<T, _deg>::degree() const {
		deg_t deg = 0;
		for (size_t i = 0; i < this->size(); i++)
			deg += (*this)[i] * (i + 1);
		return deg;
	}

	template<typename T, typename _deg>
	std::string Elementary_Symmetric_Variables<T, _deg>::name(int i, int number_of_variables) {
		return "e_" + std::to_string(i + 1);
	}


	template<typename T, typename c1, typename c2>
	Polynomial<c2> Polynomial_Basis<T, c1, c2>::operator()(Polynomial<c1> a) const {
		//set dimensions and names if nonempty
		const typename Polynomial<c2>::deg_t* gen_dims = generator_dimensions.empty() ? nullptr : generator_dimensions.data();
		const std::string* gen_names = generator_names.empty() ? nullptr : generator_names.data();
		Polynomial<c2> decomposition(number_of_variables, gen_dims, gen_names);
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

	template<typename T, typename c1, typename c2>
	Polynomial<c1> Polynomial_Basis<T, c1, c2>::operator()(const Polynomial<c2>& a) const {
		Polynomial<c1> p(number_of_variables);
		for (auto it = a.begin(); it != a.end(); ++it) {
			auto prod = compute_product(it.exponent());
			prod *= it.coeff();
			p += prod;
		}
		return p;
	}

	template<typename T, typename c1, typename c2>
	Polynomial_Basis<T, c1, c2>::Polynomial_Basis(int number_of_variables) : number_of_variables(number_of_variables) {}

	template<typename T, typename c1, typename c2>
	const auto& Polynomial_Basis<T, c1, c2>::generators() const {
		return _generators;
	}

	template<typename T, typename c1, typename c2>
	const auto& Polynomial_Basis<T, c1, c2>::dimensions() const {
		return generator_dimensions;
	}

	template<typename T, typename c1, typename c2>
	const auto& Polynomial_Basis<T, c1, c2>::names() const {
		return generator_names;
	}


	template<typename T, typename c1, typename c2>
	Polynomial<c1> Polynomial_Basis<T, c1, c2>::compute_product(const typename Polynomial<c2>::exp_t& exponent) const {
		Polynomial<c1> product(number_of_variables, 1);
		for (size_t i = 0; i < _generators.size(); i++)
			if (exponent[i] != 0)
				product *= (_generators[i] ^ exponent[i]);
		return product;
	}

	template<typename x, typename e>
	Symmetric_Basis<x, e>::Symmetric_Basis(int n) : Polynomial_Basis<Symmetric_Basis<x, e>, x, e>(n) {
		_generators.reserve(number_of_variables);
		for (int i = 1; i <= number_of_variables; i++)
			_generators.push_back(get_elementary_symmetric(i));
	}

	template<typename x, typename e>
	Polynomial<x> Symmetric_Basis<x, e>::get_elementary_symmetric(int i) const {
		Polynomial<x> poly(number_of_variables);
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

	template<typename x, typename e>
	auto Symmetric_Basis<x, e>::find_exponent(const x_t& term) const {
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
}