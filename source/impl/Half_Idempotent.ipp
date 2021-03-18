#pragma once
#include "../Half_Idempotent.hpp"

///@file
///@brief Implementation of Half_Idempotent.hpp

namespace Symmetric_Polynomials
{

	template <typename T, size_t N>
	struct Array_Vector_Wrapper : public std::array<T, N>
	{
		using std::array<T, N>::array;
		///Constructor with "size" does nothing (used to have consistent interface with vector)
		Array_Vector_Wrapper(size_t) : std::array<T, N>() {}
	};

	template <typename T>
	struct Array_Vector_Wrapper<T, 0> : public std::vector<T>
	{
		using std::vector<T>::vector;
	};

	template <typename T, typename _deg, size_t N>
	_deg Half_Idempotent_Variables<T, _deg, N>::degree() const
	{
		deg_t degree = 0;
		for (size_t i = 0; i < this->size() / 2; i++)
			degree += this->operator[](i);
		return degree;
	}

	template <typename T, typename _deg, size_t N>
	Half_Idempotent_Variables<T, _deg, N> Half_Idempotent_Variables<T, _deg, N>::operator+(const Half_Idempotent_Variables &other) const
	{
		Half_Idempotent_Variables v(this->size());
		for (size_t i = 0; i < this->size(); i++)
			v[i] = (*this)[i] + other[i];
		for (size_t i = this->size() / 2; i < this->size(); i++)
			v[i] = (v[i] > 0);
		return v;
	}

	template <typename T, typename _deg, size_t N>
	Half_Idempotent_Variables<T, _deg, N> Half_Idempotent_Variables<T, _deg, N>::operator-(const Half_Idempotent_Variables &other) const
	{
		Half_Idempotent_Variables v(this->size());
		for (size_t i = 0; i < this->size(); i++)
			v[i] = (*this)[i] - other[i];
		for (size_t i = this->size() / 2; i < this->size(); i++)
			v[i] = (v[i] > 0);
		return v;
	}

	template <typename T, typename _deg, size_t N>
	std::string Half_Idempotent_Variables<T, _deg, N>::name(int i, int number_of_variables)
	{
		if (i < number_of_variables / 2)
			return "x_" + std::to_string(i + 1);
		else
			return "y_" + std::to_string(i - number_of_variables / 2 + 1);
	}

	template <typename T, typename _deg, size_t N>
	size_t Half_Idempotent_Variables<T, _deg, N>::operator()() const
	{
		return generic_hasher(*this);
	}

	template <typename T, typename _deg>
	Twisted_Chern_Variables<T, _deg> Twisted_Chern_Variables<T, _deg>::operator+(const Twisted_Chern_Variables &other) const
	{
		Standard_Variables v(*this);
		for (size_t i = 0; i < this->size(); i++)
			v[i] += other[i];
		return v;
	}

	template <typename T, typename _deg>
	size_t Twisted_Chern_Variables<T, _deg>::operator()() const
	{
		return generic_hasher(*this);
	}

	template <typename xy, typename ch>
	Half_Idempotent_Basis<xy, ch>::Half_Idempotent_Basis(int n) : Polynomial_Basis<Half_Idempotent_Basis<xy, ch>, xy, ch>(2 * n), n(n), number_of_generators(n + (n * n + n) / 2)
	{
		set_generators();
		set_relations();
	}

	template <typename xy, typename ch>
	const auto &Half_Idempotent_Basis<xy, ch>::relations() const
	{
		return _relations;
	}

	template <typename xy, typename ch>
	const auto &Half_Idempotent_Basis<xy, ch>::generator(int s, int j) const
	{
		return _generators[index(s, j)];
	}

	template <typename xy, typename ch>
	int Half_Idempotent_Basis<xy, ch>::index(int s, int j) const
	{
		return generator_double_index.at({s, j});
	}

	template <typename xy, typename ch>
	auto Half_Idempotent_Basis<xy, ch>::create_generator(int s, int i)
	{
		Polynomial<xy> tchern(2 * n);
		xy_t mono(2 * n);
		auto c_x = Combination_Generator<typename xy_t::value_type>(n, s);
		auto c_y = Combination_Generator<typename xy_t::value_type>(n - s, i);
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
	void Half_Idempotent_Basis<xy, ch>::set_generators()
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
	void Half_Idempotent_Basis<xy, ch>::set_relations()
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
	auto Half_Idempotent_Basis<xy, ch>::find_exponent(const xy_t &term) const
	{
		chern_t exponent(number_of_generators);
		find_exponent_recursive(term, exponent);
		return exponent;
	}

	template <typename xy, typename ch>
	void Half_Idempotent_Basis<xy, ch>::find_exponent_recursive(const xy_t &term, chern_t &exponent) const
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

	///Prints all relations in the description of the fixed points of \f$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)\f$ in terms of \f$\alpha_i, c_i, \gamma_{s,j}\f$ (printed as a_i,c_i,c_{s,j} in the console)
	template <typename HIB>
	void print_half_idempotent_relations(int n, bool print, bool verify, bool verify_verbose)
	{
		HIB hib(n);
PARALLELIZE
		for (int64_t i = 0; i < (int64_t)hib.relations().size(); i++)
		{ //i is signed as MSVC still does not support OMP3.0...
			const auto &rel = hib.relations()[i];
			auto p = hib(rel);
			auto q = hib(p);
			std::stringstream ss;
			if (print)
				ss << rel << " = " << q << "\n";
			if (verify)
			{
				if (p != hib(q))
				{
					std::cerr << "Verification failed! Relation in x_i,y_i is:\n"
							  << rel << "\n Relation in a,c_i,c_{s,j} is\n"
							  << p << "\n Relation in x_i,y_i is: \n"
							  << q;
					abort();
				}
				else
				{
					ss << "Verified!";
					if (verify_verbose)
						ss << "In x, y variables both LHS and RHS are : " << p;
					ss << "\n";
				}
			}
			std::cout << ss.str();
		}
	}

}
