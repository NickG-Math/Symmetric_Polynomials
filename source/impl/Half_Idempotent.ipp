#pragma once
#include "../Half_Idempotent.hpp"

///	@file
///	@brief Implementation of Half_Idempotent.hpp

namespace symmp
{

	template <typename T, size_t N>
	struct ArrayVectorWrapper : public std::array<T, N>
	{
		using std::array<T, N>::array;
		///Constructor with "size" does nothing (used to have consistent interface with vector)
		ArrayVectorWrapper(size_t) : std::array<T, N>() {}
	};

	template <typename T>
	struct ArrayVectorWrapper<T, 0> : public std::vector<T>
	{
		using std::vector<T>::vector;
	};

	template <typename T, typename _deg, size_t N>
	_deg HalfIdempotentVariables<T, _deg, N>::degree() const
	{
		deg_t degree = 0;
		for (size_t i = 0; i < this->size() / 2; i++)
			degree += this->operator[](i);
		return degree;
	}

	template <typename T, typename _deg, size_t N>
	HalfIdempotentVariables<T, _deg, N> HalfIdempotentVariables<T, _deg, N>::operator+(const HalfIdempotentVariables &other) const
	{
		HalfIdempotentVariables v(this->size());
		for (size_t i = 0; i < this->size(); i++)
			v[i] = (*this)[i] + other[i];
		for (size_t i = this->size() / 2; i < this->size(); i++)
			v[i] = (v[i] > 0);
		return v;
	}

	template <typename T, typename _deg, size_t N>
	HalfIdempotentVariables<T, _deg, N> HalfIdempotentVariables<T, _deg, N>::operator-(const HalfIdempotentVariables &other) const
	{
		HalfIdempotentVariables v(this->size());
		for (size_t i = 0; i < this->size(); i++)
			v[i] = (*this)[i] - other[i];
		for (size_t i = this->size() / 2; i < this->size(); i++)
			v[i] = (v[i] > 0);
		return v;
	}

	template <typename T, typename _deg, size_t N>
	std::string HalfIdempotentVariables<T, _deg, N>::name(int i, int number_of_variables)
	{
		if (i < number_of_variables / 2)
			return "x_" + std::to_string(i + 1);
		else
			return "y_" + std::to_string(i - number_of_variables / 2 + 1);
	}

	template <typename T, typename _deg, size_t N>
	size_t HalfIdempotentVariables<T, _deg, N>::operator()() const
	{
		return generic_hasher(*this);
	}

	template <typename T, typename _deg>
	TwistedChernVariables<T, _deg> TwistedChernVariables<T, _deg>::operator+(const TwistedChernVariables &other) const
	{
		StandardVariables v(*this);
		for (size_t i = 0; i < this->size(); i++)
			v[i] += other[i];
		return v;
	}

	template <typename T, typename _deg>
	size_t TwistedChernVariables<T, _deg>::operator()() const
	{
		return generic_hasher(*this);
	}

	template <typename xy, typename ch>
	TwistedChernBasis<xy, ch>::TwistedChernBasis(int n) : PolynomialBasis<TwistedChernBasis<xy, ch>, xy, ch>(2 * n), n(n), number_of_generators(n + (n * n + n) / 2)
	{
		set_generators();
		set_relations();
	}

	template <typename xy, typename ch>
	const std::vector<ch> &TwistedChernBasis<xy, ch>::relations() const
	{
		return _relations;
	}

	template <typename xy, typename ch>
	const xy &TwistedChernBasis<xy, ch>::generator(int s, int j) const
	{
		return _generators[index(s, j)];
	}

	template <typename xy, typename ch>
	int TwistedChernBasis<xy, ch>::index(int s, int j) const
	{
		return generator_double_index.at({s, j});
	}

	template <typename xy, typename ch>
	auto TwistedChernBasis<xy, ch>::create_generator(int s, int i)
	{
		xy tchern;
		xy_t mono(2 * n);
		auto c_x = CombinationGenerator<typename xy_t::value_type>(n, s);
		auto c_y = CombinationGenerator<typename xy_t::value_type>(n - s, i);
		for (const auto &comb_x : c_x)
		{
			for (const auto &j : comb_x)
				mono[j] = 1;
			std::vector<typename xy_t::value_type> letters_not_in_comb_x;
			letters_not_in_comb_x.reserve(n - s);
			for (int j = 0; j < n; j++)
			{
				if (std::find(comb_x.begin(), comb_x.end(), j) == comb_x.end())
					letters_not_in_comb_x.push_back(j);
			}
			for (const auto &comb_y : c_y)
			{
				for (const auto &j : comb_y)
					mono[n + letters_not_in_comb_x[j]] = 1;
				tchern.insert(mono, 1);
				for (const auto &j : comb_y)
					mono[n + letters_not_in_comb_x[j]] = 0;
			}
			for (const auto &j : comb_x)
				mono[j] = 0;
		}
		return tchern;
	}

	template <typename xy, typename ch>
	void TwistedChernBasis<xy, ch>::set_generators()
	{
		generator_names.reserve(number_of_generators);
		generator_dimensions.reserve(number_of_generators);
		_generators.reserve(number_of_generators);
		int ind = 0;
		for (int s = 0; s <= n; s++)
		{
			for (int i = 0; i <= n - s; i++)
			{
				if (i == 0 && s == 0)
					continue;
				_generators.push_back(create_generator(s, i));
				if (i == 0)
					generator_names.push_back("c_" + std::to_string(s));
				else if (s == 0)
					generator_names.push_back("a_" + std::to_string(i));
				else
					generator_names.push_back("c_{" + std::to_string(s) + "," + std::to_string(i) + "}");
				generator_dimensions.push_back(s);
				generator_double_index[{s, i}] = ind;
				ind++;
			}
		}
	}

	template <typename xy, typename ch>
	void TwistedChernBasis<xy, ch>::set_relations()
	{
		for (int s = 0; s < n; s++)
			for (int t = s; t < n; t++)
				for (int i = 1; i <= n - s; i++)
					for (int j = 1; j <= n - t; j++)
					{
						if (t > s + i || (s == t && j < i) || (s == 0 && t != i)) //if t>s+i no relation, if s==t && j<i we have a symmetric relation, if s==0 && t!=i we can get if from the relation between a_j and a_1
							continue;
						chern_t comb(number_of_generators);
						comb[generator_double_index[{s, i}]]++;
						comb[generator_double_index[{t, j}]]++;
						_relations.emplace_back(comb, 1, generator_dimensions.data(), generator_names.data());
					}
	}

	template <typename xy, typename ch>
	auto TwistedChernBasis<xy, ch>::find_exponent(const xy_t &term) const
	{
		chern_t exponent(number_of_generators);
		find_exponent_recursive(term, exponent);
		return exponent;
	}

	template <typename xy, typename ch>
	void TwistedChernBasis<xy, ch>::find_exponent_recursive(const xy_t &term, chern_t &exponent) const
	{
		if (term[n] > 0)
		{ //clear the consecutive y_i at the start
			int consecutive_u_at_start = 0;
			xy_t us_at_start(2 * n);
			for (int i = n; i < 2 * n && term[i] > 0; i++)
			{
				us_at_start[i] = 1;
				consecutive_u_at_start++;
			}
			if (consecutive_u_at_start > 0)
				exponent[index(0, consecutive_u_at_start)] = 1;
			find_exponent_recursive(term - us_at_start, exponent);
			return;
		}
		bool found_last_u = 0;
		int last_consecutive_u_from_end = 2 * n;
		int number_consecutive_u = 0;
		for (int i = 2 * n - 1; i >= n; i--)
		{
			if (term[i] == 0 && found_last_u)
				break;
			else if (term[i] > 0)
			{
				last_consecutive_u_from_end = i;
				number_consecutive_u++;
				found_last_u = 1;
			}
		}
		if (last_consecutive_u_from_end == 2 * n)
		{ //no y_i means it's an elementary symmetric on the x_i
			for (int i = 1; i < n; i++)
				exponent[index(i, 0)] += term[i - 1] - term[i];
			exponent[index(n, 0)] += term[n - 1];
			return;
		}
		if (last_consecutive_u_from_end > n)
		{ //not first y_i so add c_{s,i} to the exponent
			auto ind = index(last_consecutive_u_from_end - n, number_consecutive_u);
			exponent[ind]++;
			find_exponent_recursive(term - _generators[ind].highest_term().exponent(), exponent);
		}
	}

	template <typename xy_poly_t, typename chern_poly_t>
	void print_half_idempotent_relations(int n, bool print, bool verify, bool verify_verbose)
	{
		TwistedChernBasis<xy_poly_t,chern_poly_t> tcb(n);
SYMMP_RUN_LOOP_IN_PARALLEL
		for (int64_t i = 0; i < (int64_t)tcb.relations().size(); i++)
		{ //i is signed as MSVC still does not support OMP3.0...
			const auto &rel = tcb.relations()[i];
			auto p = tcb(rel);
			auto q = tcb(p);
			std::stringstream ss;
			if (print)
				ss << rel << " = " << q << "\n";
			if (verify)
			{
				if (p != tcb(q))
				{
					std::cerr << "Verification failed! Relation in x_i,y_i is:\n"
							  << rel << "\n Relation in a,c_i,c_{s,j} is\n"
							  << p << "\n Relation in x_i,y_i is: \n"
							  << q;
					abort();
				}
				else
				{
					ss << "Relation verified! \n";
					if (verify_verbose)
						ss << "In x, y variables both LHS and RHS are : " << p;
					ss << "\n\n";
				}
			}
			std::cout << ss.str();
		}
	}

}
