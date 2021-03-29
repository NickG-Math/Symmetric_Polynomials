#pragma once
#include "../Polynomials.hpp"

///@file
///@brief Implementation of Polynomials.hpp

namespace symmp
{
	template <class exp_t>
	BaseContainer<exp_t>::BaseContainer(const deg_t* dimensions) : dimensions(dimensions) {}

	template <class exp_t>
	auto BaseContainer<exp_t>::compute_degree(const exp_t& exponent) const -> deg_t {
		deg_t deg;
		if constexpr (!implementation_details::has_degree_function<exp_t>::value)
		{
			if (dimensions == nullptr)
			{
				std::cerr << "If exp_t does not implement a degree() function then you must provide a valid pointer to the dimensions of the variables in the constructor";
				abort();
			}
			else
				deg = general_compute_degree<deg_t>(exponent, dimensions);
		}
		else
			deg = exponent.degree();
		return deg;
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	size_t DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::number_of_monomials() const {
		return this->size();
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	bool DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::operator==(const DefaultContainer& b) const {
		return static_cast<const data_t&>(*this) == static_cast<const data_t&>(b);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	bool DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::operator!=(const DefaultContainer& b) const {
		return !(*this == b);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::reserve(size_t n) {
		if constexpr (_ord)
			data_t::reserve(n);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::insert(const exp_t& exponent, scl_t coeff) {
		this->emplace(std::pair(this->compute_degree(exponent), exponent), coeff);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	_scl DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::ConstIterator::coeff() const {
		return data_t::const_iterator::operator*().second;
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	const _exp& DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::ConstIterator::exponent() const {
		return (data_t::const_iterator::operator*().first).second;
	}
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::ConstIterator::degree() const ->  deg_t {
		return (data_t::const_iterator::operator*().first).first;
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::ConstIterator::ConstIterator(typename data_t::const_iterator it)
		: data_t::const_iterator(it) {}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::begin() const -> ConstIterator {
		return ConstIterator(data_t::begin());
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::end() const->ConstIterator {
		return ConstIterator(data_t::end());
	}
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::highest_term() const->ConstIterator {
		if constexpr (_ord)
			return std::prev(data_t::end());
		else {
			auto it = begin();
			auto maxit = it;
			for (++it; it != end(); ++it)
				if (maxit.degree() < it.degree() || (maxit.degree() == it.degree() && maxit.exponent() < it.exponent()))
					maxit = it;
			return maxit;
		}
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	_scl& DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::Iterator::coeff() {
		return data_t::iterator::operator*().second;
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::Iterator::Iterator(typename data_t::iterator it)
		: data_t::iterator(it) {}


	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::begin() -> Iterator {
		return Iterator(data_t::begin());
	}
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	auto DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::end()->Iterator {
		return Iterator(data_t::end());
	}


	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>	
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::add(const std::pair<const std::pair<deg_t, exp_t>, scl_t>& kvp) {
		add(kvp.first, kvp.second);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::subtract(const std::pair<const std::pair<deg_t, exp_t>, scl_t>& kvp) {
		add(kvp.first, -kvp.second);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::multiply_add(const std::pair<const std::pair<deg_t, exp_t>, scl_t >& kvp1, const std::pair<const std::pair <deg_t, exp_t>, scl_t >& kvp2) {
		add(std::pair(kvp1.first.first + kvp2.first.first, kvp1.first.second + kvp2.first.second), kvp1.second * kvp2.second);
	}

	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	void DefaultContainer<_scl, _exp, _cnt, _ord, _arg...>::add(const std::pair<deg_t, exp_t> key, scl_t value) {
			auto pair = this->emplace(key, value);
			if (!pair.second)
			{ //already existing element
				pair.first->second += value;
				if (pair.first->second == 0)
					this->erase(pair.first);
			}
			//this is much faster than:	data[key] += value;	if (data[key] == 0)	data.erase(key);
		}


	template <class container_t>
	Polynomial<container_t>::Polynomial(const deg_t* dimensions, const std::string* variable_names)
		: container_t(dimensions), variable_names(variable_names) {}

	template <class container_t>
	Polynomial<container_t>::Polynomial(const exp_t& exp, scl_t coeff, const deg_t* dimensions, const std::string* variable_names)
		: Polynomial(dimensions, variable_names)
	{
		this->insert(exp, coeff);
	}

	template <class container_t>
	Polynomial<container_t>::Polynomial(int _number_of_variables, scl_t coeff, const deg_t* dimensions, const std::string* variable_names)
		: Polynomial(exp_t(_number_of_variables), coeff, dimensions, variable_names) {}


	template <class container_t>
	size_t Polynomial<container_t>::number_of_variables() const
	{
		return this->begin().exponent().size();
	}

	template <class container_t>
	auto Polynomial<container_t>::operator+=(const Polynomial& b) -> Polynomial&
	{
		this->reserve(this->number_of_monomials() + b.number_of_monomials());
		for (const auto& pair : b)
			this->add(pair);
		return *this;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator-=(const Polynomial& b) -> Polynomial&
	{
		this->reserve(this->number_of_monomials() + b.number_of_monomials());
		for (const auto& pair : b)
			this->subtract(pair);
		return *this;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator+(const Polynomial& b) const -> Polynomial
	{
		Polynomial sum;
		sum.reserve(this->number_of_monomials() + b.number_of_monomials());
		sum = *this;
		sum += b;
		return sum;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator-(const Polynomial& b) const -> Polynomial
	{
		Polynomial diff;
		diff.reserve(this->number_of_monomials() + b.number_of_monomials());
		diff = *this;
		diff -= b;
		return diff;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator*(const Polynomial& b) const -> Polynomial
	{
		Polynomial product;
		product.reserve(this->number_of_monomials() * b.number_of_monomials());
		for (const auto& paira : *this)
			for (const auto& pairb : b)
				product.multiply_add(paira, pairb);
		return product;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator*=(const Polynomial& b) -> Polynomial&
	{
		*this = operator*(b);
		return *this;
	}

	template <class container_t>
	auto Polynomial<container_t>::operator*=(scl_t coeff) -> Polynomial&
	{
		if (coeff == 0)
			*this = Polynomial(this->dimensions, variable_names);
		else if (coeff != 1)
		{
			for (auto it = this->begin(); it != this->end(); ++it)
				it.coeff() *= coeff;
		}
		return *this;
	}

	template <class container_t>
	template <typename T>
	auto Polynomial<container_t>::operator^(T p) const -> Polynomial
	{
		static_assert(std::is_integral_v<T>, "A polynomial may only be raised to a (nonnegative) integer power");
		if (p == 0)
			return Polynomial(number_of_variables(), 1, this->dimensions, variable_names);
		if (p == 1)
			return *this;
		Polynomial prod;
		// prod.reserve(power(number_of_variables(),p)); //actually slower due to how relations work
		prod = *this;
		for (T k = 1; k < p; k++)
			prod *= *this;
		return prod;
	}

	template <class container_t>
	template <typename fun>
	void Polynomial<container_t>::print(scl_t coeff, const exp_t& exponent, std::ostream& os, const fun& variable_names) const
	{
		bool havestar = 0;
		if (coeff != 1)
		{
			os << coeff;
			havestar = 1;
		}
		bool completelyzero = 1;
		for (size_t i = 0; i < exponent.size(); i++)
		{
			if (exponent[i] != 0)
			{
				completelyzero = 0;
				if (havestar)
					os << "*";
				havestar = 1;
				if (exponent[i] > 1)
					os << variable_names(i, exponent.size()) << "^" << (int)exponent[i];
				else if (exponent[i] == 1)
					os << variable_names(i, exponent.size());
			}
		}
		if (completelyzero && coeff == 1)
			os << "1";
	}

	template <class container_t>
	template <typename fun>
	void Polynomial<container_t>::print(std::ostream& os, const fun& variable_name_fun) const
	{
		auto it = this->begin();
		print(it.coeff(), it.exponent(), os, variable_name_fun);
		for (++it; it != this->end(); ++it)
		{
			os << " + ";
			print(it.coeff(), it.exponent(), os, variable_name_fun);
		}
	}

	template <class container_t>
	std::ostream& operator<<(std::ostream& os, const Polynomial<container_t>& a)
	{
		if constexpr (!implementation_details::has_name_function<typename Polynomial<container_t>::exp_t>::value)
		{
			if (a.variable_names == nullptr)
			{
				std::cerr << "If exp_t does not implement a name(int,int) function then you must provide a valid pointer to the names of the variables in the constructor";
				abort();
			}
			else
				a.print(os, [=](int i, int) { return a.variable_names[i]; });
		}
		else
			a.print(os, Polynomial<container_t>::exp_t::name);
		return os;
	}
}
