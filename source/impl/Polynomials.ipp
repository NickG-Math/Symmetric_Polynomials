#pragma once
#include "../Polynomials.hpp"

///@file
///@brief Implementation of Polynomials.hpp

namespace
{ //SFINAE stuff

	template <typename T>
	static constexpr std::false_type test_degree_existence(...);

	template <typename T>
	static constexpr decltype(std::declval<T>().degree(), std::true_type()) test_degree_existence(int);

	template <typename T>
	using has_degree_function = decltype(test_degree_existence<T>(0));

	template <typename T>
	static constexpr std::false_type test_name_existence(...);

	template <typename T>
	static constexpr decltype(std::declval<T>().name(std::declval<int>(), std::declval<int>()), std::true_type()) test_name_existence(int);

	template <typename T>
	using has_name_function = decltype(test_name_existence<T>(0));

	template <typename T>
	static constexpr std::false_type test_reserve_existence(...);

	template <typename T>
	static constexpr decltype(std::declval<T>().reserve(std::declval<size_t>()), std::true_type()) test_reserve_existence(int);

	template <typename T>
	using has_reserve_function = decltype(test_reserve_existence<T>(0));

	///Pairwise addition of pairs
	template <typename T, typename S>
	std::pair<T, S> operator+(const std::pair<T, S> &a, const std::pair<T, S> b)
	{
		return std::pair<T, S>(a.first + b.first, a.second + b.second);
	}
	///Given a pair, hash only the second parameter
	template <typename T, typename S>
	struct hash_only_second_in_pair
	{
		///Hash only second parameter
		auto operator()(const std::pair<T, S> &pair) const
		{
			return pair.second();
		}
	};

	// //naive power, should be fine
	// size_t power(size_t a, size_t b){
	// 	size_t p=1;
	// 	for (size_t i=0; i<b; i++)
	// 		p*=a;
	// 	return p;
	// }
}

namespace Symmetric_Polynomials
{

	///Ordered specialization of default container
	template <typename scl_t, typename exp_t>
	struct default_container<scl_t, exp_t, 1> : public std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t>
	{
		using std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t>::map;
		static constexpr bool ordered = 1; ///<This is needed to tell the polynomial that the container is ordered
	};

	///Unordered specialization of default container
	template <typename scl_t, typename exp_t>
	struct default_container<scl_t, exp_t, 0> : public std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hash_only_second_in_pair<typename exp_t::deg_t, exp_t>>
	{
		using std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hash_only_second_in_pair<typename exp_t::deg_t, exp_t>>::unordered_map;
		static constexpr bool ordered = 0; ///<This is needed to tell the polynomial that the container is unordered
	};

	template <typename container_t>
	Polynomial<container_t>::Polynomial(int _number_of_variables, const deg_t *_dimensions, const std::string *_variable_names)
		: _number_of_variables(_number_of_variables), _dimensions(_dimensions), _variable_names(_variable_names) {}

	template <typename container_t>
	Polynomial<container_t>::Polynomial(const exp_t &exp, scl_t coeff, const deg_t *_dimensions, const std::string *_variable_names)
		: Polynomial(exp.size(), _dimensions, _variable_names)
	{
		insert(exp, coeff);
	}

	template <typename container_t>
	Polynomial<container_t>::Polynomial(int _number_of_variables, scl_t coeff, const deg_t *_dimensions, const std::string *_variable_names)
		: Polynomial(exp_t(_number_of_variables), coeff, _dimensions, _variable_names) {}

	template <typename container_t>
	auto Polynomial<container_t>::number_of_variables() const
	{
		return _number_of_variables;
	}

	template <typename container_t>
	void Polynomial<container_t>::reserve(size_t n)
	{
		if constexpr (has_reserve_function<container_t>::value)
			data.reserve(n);
	}

	template <typename container_t>
	auto Polynomial<container_t>::number_of_monomials() const
	{
		return data.size();
	}

	template <typename container_t>
	void Polynomial<container_t>::clear()
	{
		data.clear();
	}

	template <typename container_t>
	void Polynomial<container_t>::insert(const exp_t &exponent, scl_t scalar)
	{
		data.emplace(typename container_t::key_type(compute_degree(exponent), exponent), scalar);
	}

	template <typename container_t>
	auto Polynomial<container_t>::const_iterator::coeff() const
	{
		return it->second;
	}

	template <typename container_t>
	const auto &Polynomial<container_t>::const_iterator::exponent() const
	{
		return (it->first).second;
	}

	template <typename container_t>
	auto Polynomial<container_t>::const_iterator::degree() const
	{
		return (it->first).first;
	}

	template <typename container_t>
	typename Polynomial<container_t>::const_iterator &Polynomial<container_t>::const_iterator::operator++()
	{
		++it;
		return *this;
	}

	template <typename container_t>
	bool Polynomial<container_t>::const_iterator::operator==(const_iterator second) const
	{
		return (it == second.it);
	}

	template <typename container_t>
	bool Polynomial<container_t>::const_iterator::operator!=(const_iterator second) const
	{
		return (it != second.it);
	}

	template <typename container_t>
	Polynomial<container_t>::const_iterator::const_iterator(typename container_t::const_iterator it) : it(it) {}

	template <typename container_t>
	typename Polynomial<container_t>::const_iterator Polynomial<container_t>::begin() const
	{
		return const_iterator(data.begin());
	}

	template <typename container_t>
	typename Polynomial<container_t>::const_iterator Polynomial<container_t>::end() const
	{
		return const_iterator(data.end());
	}

	template <typename container_t>
	typename Polynomial<container_t>::const_iterator Polynomial<container_t>::highest_term() const
	{
		if constexpr (container_t::ordered)
			return std::prev(data.end());
		else
		{
			auto it = begin();
			auto maxit = it;
			for (++it; it != end(); ++it)
			{
				if (maxit.degree() < it.degree() || (maxit.degree() == it.degree() && maxit.exponent() < it.exponent()))
					maxit = it;
			}
			return maxit;
		}
	}

	template <typename container_t>
	Polynomial<container_t> &Polynomial<container_t>::operator+=(const Polynomial &b)
	{
		reserve(number_of_monomials() + b.number_of_monomials());
		for (const auto &pair : b.data)
			insert_add_erase(pair.first, pair.second);
		return *this;
	}

	template <typename container_t>
	Polynomial<container_t> &Polynomial<container_t>::operator-=(const Polynomial &b)
	{
		reserve(number_of_monomials() + b.number_of_monomials());
		for (const auto &pair : b.data)
			insert_add_erase(pair.first, -pair.second);
		return *this;
	}

	template <typename container_t>
	Polynomial<container_t> Polynomial<container_t>::operator+(const Polynomial &b) const
	{
		Polynomial sum;
		sum.reserve(number_of_monomials() + b.number_of_monomials());
		sum = *this;
		sum += b;
		return sum;
	}

	template <typename container_t>
	Polynomial<container_t> Polynomial<container_t>::operator-(const Polynomial &b) const
	{
		Polynomial diff;
		diff.reserve(number_of_monomials() + b.number_of_monomials());
		diff = *this;
		diff -= b;
		return diff;
	}

	template <typename container_t>
	Polynomial<container_t> Polynomial<container_t>::operator*(const Polynomial &b) const
	{
		Polynomial product(_number_of_variables);
		product.reserve(number_of_monomials() * b.number_of_monomials());
		for (const auto &paira : data)
			for (const auto &pairb : b.data)
				product.insert_add_erase(paira.first + pairb.first, paira.second * pairb.second);
		return product;
	}

	template <typename container_t>
	Polynomial<container_t> &Polynomial<container_t>::operator*=(const Polynomial &b)
	{
		*this = operator*(b);
		return *this;
	}

	template <typename container_t>
	Polynomial<container_t> &Polynomial<container_t>::operator*=(scl_t coeff)
	{
		if (coeff == 0)
			clear();
		else if (coeff != 1)
		{
			for (const auto &paira : data)
				data[paira.first] *= coeff;
		}
		return *this;
	}

	template <typename container_t>
	template <typename any_int_type>
	Polynomial<container_t> Polynomial<container_t>::operator^(any_int_type p) const
	{
		static_assert(std::is_integral_v<any_int_type>, "A polynomial may only be raised to a (nonnegative) integer power");
		if (p == 0)
			return Polynomial(_number_of_variables, 1, _dimensions, _variable_names);
		if (p == 1)
			return *this;
		Polynomial prod;
		// prod.reserve(power(number_of_variables(),p)); //actually slower due to how relations work
		prod = *this;
		for (any_int_type k = 1; k < p; k++)
			prod *= *this;
		return prod;
	}

	template <typename container_t>
	std::string Polynomial<container_t>::print() const
	{
		if constexpr (!has_name_function<exp_t>::value)
		{
			if (_variable_names == nullptr)
			{
				std::cerr << "If the exp_t class does not implement a degree() function then you must provide a valid pointer to the dimensions of the variables in the polynomial constructor";
				abort();
			}
			else
				return print([=](int i, int) { return _variable_names[i]; });
		}
		else
			return print(exp_t::name);
	}

	template <typename container_t>
	bool Polynomial<container_t>::operator==(const Polynomial &b) const
	{
		return (data == b.data);
	}

	template <typename container_t>
	bool Polynomial<container_t>::operator!=(const Polynomial &b) const
	{
		return !(*this == b);
	}

	template <typename container_t>
	auto Polynomial<container_t>::compute_degree(const exp_t &exponent)
	{
		deg_t deg;
		if constexpr (!has_degree_function<exp_t>::value)
		{
			if (_dimensions == nullptr)
			{
				std::cerr << "If the exp_t class does not implement a degree() function then you must provide a valid pointer to the dimensions of the variables in the polynomial constructor";
				abort();
			}
			else
				deg = general_compute_degree<deg_t>(exponent, _dimensions);
		}
		else
			deg = exponent.degree();
		return deg;
	}

	template <typename container_t>
	template <typename key_t>
	void Polynomial<container_t>::insert_add_erase(const key_t &key, scl_t value)
	{
		auto pair = data.emplace(key, value);
		if (!pair.second)
		{ //already existing element
			pair.first->second += value;
			if (pair.first->second == 0)
				data.erase(pair.first);
		}
		//this is much faster than:	data[key] += value;	if (data[key] == 0)	data.erase(key);
	}

	template <typename container_t>
	template <typename fun>
	void Polynomial<container_t>::print(scl_t coeff, const exp_t &exponent, std::stringstream &ss, const fun &_variable_names) const
	{
		bool havestar = 0;
		if (coeff != 1)
		{
			ss << coeff;
			havestar = 1;
		}
		bool completelyzero = 1;
		for (size_t i = 0; i < exponent.size(); i++)
		{
			if (exponent[i] != 0)
			{
				completelyzero = 0;
				if (havestar)
					ss << "*";
				havestar = 1;
				if (exponent[i] > 1)
					ss << _variable_names(i, _number_of_variables) << "^" << (int)exponent[i];
				else if (exponent[i] == 1)
					ss << _variable_names(i, _number_of_variables);
			}
		}
		if (completelyzero && coeff == 1)
			ss << "1";
	}

	template <typename container_t>
	template <typename fun>
	std::string Polynomial<container_t>::print(const fun &variable_name_fun) const
	{
		std::stringstream ss;
		auto it = begin();
		print(it.coeff(), it.exponent(), ss, variable_name_fun);
		for (++it; it != end(); ++it)
		{
			ss << " + ";
			print(it.coeff(), it.exponent(), ss, variable_name_fun);
		}
		return ss.str();
	}

	template <typename container_t>
	std::ostream &operator<<(std::ostream &os, const Polynomial<container_t> &a)
	{
		os << a.print();
		return os;
	}

}
