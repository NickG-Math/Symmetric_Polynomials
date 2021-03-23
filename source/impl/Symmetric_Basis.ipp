#pragma once
#include "../Symmetric_Basis.hpp"

///@file
///@brief Implementation of SymmetricBasis.hpp

namespace symmp
{

	template <typename T, typename _deg>
	_deg StandardVariables<T, _deg>::degree() const
	{
		deg_t sum = 0;
		for (const auto i : *this)
			sum += i;
		return sum;
	}

	template <typename T, typename _deg>
	std::string StandardVariables<T, _deg>::name(int i, int)
	{
		return "x_" + std::to_string(i + 1);
	}

	template <typename T, typename _deg>
	StandardVariables<T, _deg> StandardVariables<T, _deg>::operator+(const StandardVariables &other) const
	{
		StandardVariables v(*this);
		for (size_t i = 0; i < this->size(); i++)
			v[i] += other[i];
		return v;
	}

	template <typename T, typename _deg>
	size_t StandardVariables<T, _deg>::operator()() const
	{
		return generic_hasher(*this);
	}

	template <typename T, typename _deg>
	_deg ElementarySymmetricVariables<T, _deg>::degree() const
	{
		deg_t deg = 0;
		for (size_t i = 0; i < this->size(); i++)
			deg += (*this)[i] * (i + 1);
		return deg;
	}

	template <typename T, typename _deg>
	std::string ElementarySymmetricVariables<T, _deg>::name(int i, int)
	{
		return "e_" + std::to_string(i + 1);
	}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	new_poly_t PolynomialBasis<T, orig_poly_t, new_poly_t>::operator()(orig_poly_t a) const
	{
		//set dimensions and names if nonempty
		const typename new_poly_t::deg_t *gen_dims = generator_dimensions.empty() ? nullptr : generator_dimensions.data();
		const std::string *gen_names = generator_names.empty() ? nullptr : generator_names.data();
		new_poly_t decomposition(gen_dims, gen_names);
		while (true)
		{
			auto max = a.highest_term();
			auto exponent = static_cast<const T *>(this)->find_exponent(max.exponent());
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

	template <typename T, typename orig_poly_t, typename new_poly_t>
	orig_poly_t PolynomialBasis<T, orig_poly_t, new_poly_t>::operator()(const new_poly_t&a) const
	{
		orig_poly_t p;
		for (auto it = a.begin(); it != a.end(); ++it)
		{
			auto prod = compute_product(it.exponent());
			prod *= it.coeff();
			p += prod;
		}
		return p;
	}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	PolynomialBasis<T, orig_poly_t, new_poly_t>::PolynomialBasis(int number_of_variables) : number_of_variables(number_of_variables) {}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	const std::vector<orig_poly_t> & PolynomialBasis<T, orig_poly_t, new_poly_t>::generators() const
	{
		return _generators;
	}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	const std::vector<typename new_poly_t::deg_t> & PolynomialBasis<T, orig_poly_t, new_poly_t>::dimensions() const
	{
		return generator_dimensions;
	}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	const std::vector<std::string>& PolynomialBasis<T, orig_poly_t, new_poly_t>::names() const
	{
		return generator_names;
	}

	template <typename T, typename orig_poly_t, typename new_poly_t>
	orig_poly_t PolynomialBasis<T, orig_poly_t, new_poly_t>::compute_product(const typename new_poly_t::exp_t &exponent) const
	{
		orig_poly_t product(number_of_variables, 1);
		for (size_t i = 0; i < _generators.size(); i++)
			if (exponent[i] != 0)
				product *= (_generators[i] ^ exponent[i]);
		return product;
	}

	template <typename x_poly_t, typename e_poly_t>
	SymmetricBasis<x_poly_t, e_poly_t>::SymmetricBasis(int n) : PolynomialBasis<SymmetricBasis<x_poly_t, e_poly_t>, x_poly_t, e_poly_t>(n)
	{
		_generators.reserve(number_of_variables);
		for (int i = 1; i <= number_of_variables; i++)
			_generators.push_back(get_elementary_symmetric(i));
	}

	template <typename x_poly_t, typename e_poly_t>
	x_poly_t SymmetricBasis<x_poly_t, e_poly_t>::get_elementary_symmetric(int i) const
	{
		x_poly_t poly;
		x_t mono(number_of_variables);
		auto c = CombinationGenerator(number_of_variables, i);
		for (const auto &comb : c)
		{
			for (const auto j : comb)
				mono[j] = 1;
			poly.insert(mono, 1);
			for (const auto j : comb)
				mono[j] = 0;
		}
		return poly;
	}

	template <typename x_poly_t, typename e_poly_t>
	auto SymmetricBasis<x_poly_t, e_poly_t>::find_exponent(const x_t &term) const -> e_t
	{
		e_t exponent;
		exponent.reserve(number_of_variables);
		for (size_t i = 0; i < term.size(); i++)
		{
			if (i+1 == term.size())
				exponent.push_back(term[i]);
			else
				exponent.push_back(term[i] - term[i + 1]);
		}
		return exponent;
	}
}