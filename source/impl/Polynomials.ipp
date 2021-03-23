#pragma once
#include "../Polynomials.hpp"

///@file
///@brief Implementation of Polynomials.hpp

namespace symmp
{
	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	Polynomial<t1, t2, t3, t4, t5...>::Polynomial(const deg_t* _dimensions, const std::string* _variable_names)
		: _dimensions(_dimensions), _variable_names(_variable_names) {}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	Polynomial<t1, t2, t3, t4, t5...>::Polynomial(const exp_t& exp, scl_t coeff, const deg_t* _dimensions, const std::string* _variable_names)
		: Polynomial(_dimensions, _variable_names)
	{
		insert(exp, coeff);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	Polynomial<t1, t2, t3, t4, t5...>::Polynomial(int _number_of_variables, scl_t coeff, const deg_t* _dimensions, const std::string* _variable_names)
		: Polynomial(exp_t(_number_of_variables), coeff, _dimensions, _variable_names) {}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	void Polynomial<t1, t2, t3, t4, t5...>::reserve(size_t n)
	{
		if constexpr (implementation_details::has_reserve_function<data_t>::value)
			data.reserve(n);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	size_t Polynomial<t1, t2, t3, t4, t5...>::number_of_variables() const
	{
		return begin().exponent().size();
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	size_t Polynomial<t1, t2, t3, t4, t5...>::number_of_monomials() const
	{
		return data.size();
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	void Polynomial<t1, t2, t3, t4, t5...>::insert(const exp_t& exponent, scl_t scalar)
	{
		data.emplace(std::pair(compute_degree(exponent), exponent), scalar);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::constIterator::coeff() const -> scl_t
	{
		return it->second;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::constIterator::exponent() const -> const exp_t&
	{
		return (it->first).second;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::constIterator::degree() const -> deg_t
	{
		return (it->first).first;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::constIterator::operator++() -> constIterator&
	{
		++it;
		return *this;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	bool Polynomial<t1, t2, t3, t4, t5...>::constIterator::operator==(constIterator second) const
	{
		return (it == second.it);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	bool Polynomial<t1, t2, t3, t4, t5...>::constIterator::operator!=(constIterator second) const
	{
		return (it != second.it);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	Polynomial<t1, t2, t3, t4, t5...>::constIterator::constIterator(it_t it) : it(it) {}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::begin() const ->constIterator
	{
		return constIterator(data.begin());
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::end() const ->constIterator
	{
		return constIterator(data.end());
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::highest_term() const ->constIterator
	{
		if constexpr (t4) //if container is ordered
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

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator+=(const Polynomial& b) -> Polynomial&
	{
		reserve(number_of_monomials() + b.number_of_monomials());
		for (const auto& pair : b.data)
			insert_add_erase(pair.first, pair.second);
		return *this;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator-=(const Polynomial& b) -> Polynomial&
	{
		reserve(number_of_monomials() + b.number_of_monomials());
		for (const auto& pair : b.data)
			insert_add_erase(pair.first, -pair.second);
		return *this;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator+(const Polynomial& b) const -> Polynomial
	{
		Polynomial sum;
		sum.reserve(number_of_monomials() + b.number_of_monomials());
		sum = *this;
		sum += b;
		return sum;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator-(const Polynomial& b) const -> Polynomial
	{
		Polynomial diff;
		diff.reserve(number_of_monomials() + b.number_of_monomials());
		diff = *this;
		diff -= b;
		return diff;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator*(const Polynomial& b) const -> Polynomial
	{
		Polynomial product;
		product.reserve(number_of_monomials() * b.number_of_monomials());
		for (const auto& paira : data)
			for (const auto& pairb : b.data)
				product.insert_add_erase(implementation_details::operator+(paira.first, pairb.first), paira.second * pairb.second);
		return product;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator*=(const Polynomial& b) -> Polynomial&
	{
		*this = operator*(b);
		return *this;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator*=(scl_t coeff) -> Polynomial&
	{
		if (coeff == 0)
			*this=Polynomial(_dimensions,_variable_names);
		else if (coeff != 1)
		{
			for (const auto& paira : data)
				data[paira.first] *= coeff;
		}
		return *this;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	template <typename T>
	auto Polynomial<t1, t2, t3, t4, t5...>::operator^(T p) const -> Polynomial
	{
		static_assert(std::is_integral_v<T>, "A polynomial may only be raised to a (nonnegative) integer power");
		if (p == 0)
			return Polynomial(number_of_variables(), 1, _dimensions, _variable_names);
		if (p == 1)
			return *this;
		Polynomial prod;
		// prod.reserve(power(number_of_variables(),p)); //actually slower due to how relations work
		prod = *this;
		for (T k = 1; k < p; k++)
			prod *= *this;
		return prod;
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	std::string Polynomial<t1, t2, t3, t4, t5...>::print() const
	{
		if constexpr (!implementation_details::has_name_function<exp_t>::value)
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

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	bool Polynomial<t1, t2, t3, t4, t5...>::operator==(const Polynomial& b) const
	{
		return (data == b.data);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	bool Polynomial<t1, t2, t3, t4, t5...>::operator!=(const Polynomial& b) const
	{
		return !(*this == b);
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	auto Polynomial<t1, t2, t3, t4, t5...>::compute_degree(const exp_t& exponent)->deg_t
	{
		deg_t deg;
		if constexpr (!implementation_details::has_degree_function<exp_t>::value)
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

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	template <typename key_t>
	void Polynomial<t1, t2, t3, t4, t5...>::insert_add_erase(const key_t& key, scl_t value)
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

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	template <typename fun>
	void Polynomial<t1, t2, t3, t4, t5...>::print(scl_t coeff, const exp_t& exponent, std::stringstream& ss, const fun& _variable_names) const
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
					ss << _variable_names(i, exponent.size()) << "^" << (int)exponent[i];
				else if (exponent[i] == 1)
					ss << _variable_names(i, exponent.size());
			}
		}
		if (completelyzero && coeff == 1)
			ss << "1";
	}

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	template <typename fun>
	std::string Polynomial<t1, t2, t3, t4, t5...>::print(const fun& variable_name_fun) const
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

	template <typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
	std::ostream& operator<<(std::ostream& os, const Polynomial<t1, t2, t3, t4, t5...>& a)
	{
		os << a.print();
		return os;
	}

}
