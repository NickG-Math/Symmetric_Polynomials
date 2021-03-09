#pragma once
#include <functional>
#include <map>
#include <unordered_map>
#include "General.h"

///@file
///@brief Contains the class of polynomials in multiple variables. 


namespace { //SFINAE stuff

	template<typename T>
	static constexpr std::false_type test_degree_existence(...);

	template<typename T>
	static constexpr decltype(std::declval<T>().degree(), std::true_type()) test_degree_existence(int);

	template<typename T>
	using has_degree_function = decltype(test_degree_existence<T>(0));

	template<typename T>
	static constexpr std::false_type test_name_existence(...);

	template<typename T>
	static constexpr decltype(std::declval<T>().name(std::declval<int>(), std::declval<int>()), std::true_type()) test_name_existence(int);

	template<typename T>
	using has_name_function = decltype(test_name_existence<T>(0));

	template<typename T>
	static constexpr std::false_type test_reserve_existence(...);

	template<typename T>
	static constexpr decltype(std::declval<T>().reserve(std::declval<size_t>()), std::true_type()) test_reserve_existence(int);

	template<typename T>
	using has_reserve_function = decltype(test_reserve_existence<T>(0));

}

namespace Symmetric_Polynomials { 

	///Pairwise addition of pairs
	template<typename T, typename S>
	std::pair<T, S> operator+(const std::pair<T, S>& a, const std::pair<T, S> b) {
		return std::pair<T, S>(a.first + b.first, a.second + b.second);
	}

	///Given a pair, hash only the second parameter
	template<typename T, typename S>
	struct hash_only_second_in_pair {
		///Hash only second parameter
		auto operator()(const std::pair<T, S>& pair) const {
			return pair.second();
		}
	};

	////////////////////////////////////////////////////////////////////////
	///The default container for polynomials given scalar (coefficient) type (eg int,float), exponent type (eg Standard_Variables , Half_Idempotent_Variables) and whether the container should be ordered or not
	//
	///When ordered=1 this is a wrapper of std::map and when ordered=0 this is a wrapper of std::unordered_map
	////////////////////////////////////////////////////////////////////////
	template<typename scl_t, typename exp_t, bool ordered = 1>
	struct default_container;

	///Ordered specialization of default container
	template<typename scl_t, typename exp_t>
	struct default_container<scl_t, exp_t, 1> : public std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t> {
		using std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t>::map;
		static constexpr bool ordered = 1; ///<This is needed to tell the polynomial that the container is ordered
	};

	///Unordered specialization of default container
	template<typename scl_t, typename exp_t>
	struct default_container<scl_t, exp_t, 0> : public std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hash_only_second_in_pair<typename exp_t::deg_t, exp_t>> {
		using std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hash_only_second_in_pair<typename exp_t::deg_t, exp_t>>::unordered_map;
		static constexpr bool ordered = 0;  ///<This is needed to tell the polynomial that the container is unordered
	};

	///////////////////////////////////////////////////////////////////////////////////
	///Class of Polynomials in multiple variables
	//
	///The template parameter container_t specifies the container for the monomials in the polynomial eg default_container<int, Standard_Variables<>, 1> or default_container<float, Half_Idempotent_Variables<>, 0>
	///
	///
	///In essence the container must map pairs (degree,exponent)=\f$(d,e)\f$ to scalars=\f$c\f$ representing the monomial \f$cx_1^{e_1}\cdots x_n^{e_n}\f$ of degree \f$d\f$.
	///
	///
	///Requirements of container_t:
	///
	///- container_t must have a static constexpr bool variable ordered that specifies whether it keeps the monomials in increasing order (yes if ordered=1, no if ordered=0)
	///
	///- If ordered=1, container_t should be equivalent to std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t> 
	///
	///- If ordered=0, container_t should be equivalent to std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hasher>
	///
	///In particular container_t must also have: operator[](), const_iterator with begin(),end() methods, and methods insert() and erase().
	///
	///Finally, it must have the typedefs mapped_type (the scalar), key_type::first_type (the degree), key_type::second_type (the exponent)
	//
	///The exponent type key_type::second_type, henceforth exp_t, must have functionality similar to Standard_Variables or Half_Idempotent_Variables
	///
	///
	///Requirements of exp_t (i.e. key_type::second_type):
	///
	///- A typedef deg_t that represents the degree type (eg int or uint64_t)
	///
	///- Basic vector functionality (constructor that takes int n and produces exponent of 0's with that number of variables n, operator[])
	///
	///- It should have an operator + to be used in the product of monomials and an operator - to be used in the division of monomials (if products/divisions need to be used).
	///
	///- If ordered=1 it needs to have an operator < (implementing a lexicographic ordering) and if ordered=0 it needs to have an operator () implementing a hash function
	///
	///- It may also optionally have a degree() function and a static name(int,int) function that compute the degree of the exponent and print the names of the variables (first parameter is the index of the variable, second is the total number of variables).
	///
	///- if exp_t does not have a degree() function then the user needs to provide the dimensions of the variables through a deg_t* (the user also needs to make sure that this pointer is not invalidated as long as the polynomial is used)
	///
	///- if exp_t does not have a static name(int,int) function then the user needs to provide the names of the variables through an std::string* (the user also user needs to make sure that this pointer is not invalidated as long as the polynomial is used)
	///////////////////////////////////////////////////////////////////////////////////
	template<typename container_t>
	class Polynomial {
	public:

		typedef typename container_t::mapped_type scl_t; ///<The scalar type
		typedef typename container_t::key_type::second_type exp_t; ///<The exponent type
		typedef typename container_t::key_type::first_type deg_t; ///<The degree type

		///Default constructor (empty polynomial)
		Polynomial() {}

		///Constructs empty polynomial using number of variables and optionally dimensions of the variables and names of the variables (used only when exp_t does not implement degree() and/or name(int,int))
		Polynomial(int _number_of_variables, const deg_t* _dimensions = nullptr, const std::string* _variable_names = nullptr) : _number_of_variables(_number_of_variables), _dimensions(_dimensions), _variable_names(_variable_names) {}

		///Constructs polynomial using monomial (in the form of exponent + coefficient) and optionally dimensions of the variables and names of the variables provided as pointers (used only when exp_t does not implement degree() and/or name(int,int))
		Polynomial(const exp_t& exp, scl_t coeff, const deg_t* _dimensions = nullptr, const std::string* _variable_names = nullptr) : Polynomial(exp.size(), _dimensions, _variable_names) {
			insert(exp, coeff);
		}

		///Constructs constant polynomial using coefficient, number of variables and optionally dimensions of the variables and names of the variables (used only when exp_t does not implement degree() and/or name(int,int))
		Polynomial(int _number_of_variables, scl_t coeff, const deg_t* _dimensions = nullptr, const std::string* _variable_names = nullptr) : Polynomial(exp_t(_number_of_variables), coeff, _dimensions, _variable_names) {}

		///Returns the numbers of variables of the polynomial
		auto number_of_variables() const { return _number_of_variables; }

		///Reserve number of monomials (does nothing for ordered polynomials)
		void reserve(size_t n) {
			if constexpr (has_reserve_function<container_t>::value)
				data.reserve(n);
		}

		///Returns number of monomials in polynomial
		auto number_of_monomials() const { return data.size(); }

		///Clears all data of polynomial apart from number of variables
		void clear() { data.clear(); }

		///Inserts monomial given exponent and coefficient (it is the user's responsibility to make sure that the coefficient is nonzero). If a monomial with same exponent already exists then this does nothing
		void insert(const exp_t& exponent, scl_t scalar) {
			data.emplace(typename container_t::key_type(compute_degree(exponent), exponent), scalar);
		}

		///Constant iterator going through the monomials of the polynomial in increasing order
		class const_iterator {
		public:
			///Returns the coefficient of the monomial
			auto coeff() const {
				return it->second;
			}
			///Returns the exponent of the monomial
			const exp_t& exponent() const {
				return (it->first).second;
			}
			///Returns the degree of the monomial
			auto degree() const {
				return (it->first).first;
			}
			///Increases iterator
			const_iterator& operator ++() {
				++it;
				return *this;
			}
			///Equality of iterators
			bool operator ==(const_iterator second) const {
				return (it == second.it);
			}
			///Inequality of iterators
			bool operator !=(const_iterator second) const {
				return (it != second.it);
			}
		private:
			typename container_t::const_iterator it;
			const_iterator(typename container_t::const_iterator it) : it(it) {}
			friend class Polynomial;
		};

		///Starting iterator
		const_iterator begin() const {
			return const_iterator(data.begin());
		}
		///Ending iterator
		const_iterator end() const {
			return const_iterator(data.end());
		}

		///Returns highest term monomial in a polynomial
		const_iterator highest_term() const {
			if constexpr (container_t::ordered)
				return std::prev(data.end());
			else {
				auto it = begin();
				auto maxit = it;
				for (++it; it != end(); ++it) {
					if (maxit.degree() < it.degree() || (maxit.degree() == it.degree() && maxit.exponent() < it.exponent()))
						maxit = it;
				}
				return maxit;
			}

		}

		///Adds polynomial on the RHS to polynomial on the LHS (in place)
		Polynomial& operator+=(const Polynomial& b) {
			reserve(number_of_monomials() + b.number_of_monomials());
			for (const auto& pair : b.data)
				insert_add_erase(pair.first, pair.second);
			return *this;
		}

		///Subtracts polynomial on the RHS from polynomial on the LHS (in place)
		Polynomial& operator-=(const Polynomial& b) {
			reserve(number_of_monomials() + b.number_of_monomials());
			for (const auto& pair : b.data)
				insert_add_erase(pair.first, -pair.second);
			return *this;
		}

		///Adds polynomials
		Polynomial operator+(const Polynomial& b) const {
			Polynomial sum;
			sum.reserve(number_of_monomials() + b.number_of_monomials());
			sum = *this;
			sum += b;
			return sum;
		}

		///Subtracts polynomials
		Polynomial operator-(const Polynomial& b) const {
			Polynomial diff;
			diff.reserve(number_of_monomials() + b.number_of_monomials());
			diff = *this;
			diff -= b;
			return diff;
		}

		///Multiplies polynomials
		Polynomial operator* (const Polynomial& b) const {
			Polynomial product(_number_of_variables);
			product.reserve(number_of_monomials() * b.number_of_monomials());
			for (const auto& paira : data)
				for (const auto& pairb : b.data)
					product.insert_add_erase(paira.first + pairb.first, paira.second * pairb.second);
			return product;
		}

		///Multiplies polynomial on the LHS with polynomial on the RHS (not in place)
		Polynomial& operator*= (const Polynomial& b) {
			*this = operator*(b);
			return *this;
		}

		///Multiplies polynomial on the LHS with scalar coefficient on the RHS (in place)
		Polynomial& operator*= (scl_t coeff) {
			if (coeff == 0)
				clear();
			else if (coeff != 1) {
				for (const auto& paira : data)
					data[paira.first] *= coeff;
			}
			return *this;
		}

		///Raises polynomial to given power (does nothing if p<0)
		template<typename any_int_type=int>
		Polynomial operator^(any_int_type p) const {
			static_assert(std::is_integral_v<any_int_type>, "A polynomial may only be raised to a (nonnegative) integer power");
			if (p == 0)
				return Polynomial(_number_of_variables, 1, _dimensions, _variable_names);
			if (p == 1)
				return *this;
			auto prod = *this;
			for (any_int_type k = 1; k < p; k++)
				prod *= *this;
			return prod;
		}

		///Prints polynomial
		std::string print() const {
			if constexpr (!has_name_function<exp_t>::value) {
				if (_variable_names == nullptr)
					throw("If the exp_t class does not implement a name(int,int) function then you must provide a valid pointer to the names of the variables in the polynomial constructor");
				else
					return print([=](int i, int no) {return _variable_names[i];});
			}
			else
				return print(exp_t::name);
		}

		///Equality of polynomials
		bool operator==(const Polynomial& b) const {
			return (data == b.data);
		}

		///Inequality of polynomials
		bool operator!=(const Polynomial& b)  const {
			return !(*this == b);
		}

	private:
	
		int _number_of_variables;
		container_t data;
		const deg_t* _dimensions; //if exp_t does not have the appropriate method
		const std::string* _variable_names; //if exp_t does not have the appropriate method




		///Finds degree of given exponent
		deg_t compute_degree(const exp_t& exponent) {
			deg_t deg;
			if constexpr (!has_degree_function<exp_t>::value) {
				if (_dimensions == nullptr)
					throw("If the exp_t class does not implement a degree() function then you must provide a valid pointer to the dimensions of the variables in the polynomial constructor");
				else
					deg = general_compute_degree<deg_t>(exponent, _dimensions);
			}
			else
				deg = exponent.degree();
			return deg;
		}

		template<typename key_t>
		void insert_add_erase(const key_t& key, scl_t value) {
			auto pair = data.emplace(key, value);
			if (!pair.second) {//already existing element
				pair.first->second += value;
				if (pair.first->second == 0)
					data.erase(pair.first);
			}
			//this is much faster than:	data[key] += value;	if (data[key] == 0)	data.erase(key);
		}

		///Print monomial using given variable names
		void print(scl_t coeff, const exp_t& exponent, std::stringstream& ss, const std::function<std::string(int, int)>& _variable_names) const {
			bool havestar = 0;
			if (coeff != 1) {
				ss << coeff;
				havestar = 1;
			}
			bool completelyzero = 1;
			for (size_t i = 0; i < exponent.size(); i++) {
				if (exponent[i] != 0) {
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


		///Print polynomial using given variable names
		std::string print(const std::function<std::string(int, int)>& variable_name_fun) const {
			std::stringstream ss;
			auto it = begin();
			print(it.coeff(), it.exponent(), ss, variable_name_fun);
			for (++it; it != end(); ++it) {
				ss << " + ";
				print(it.coeff(), it.exponent(), ss, variable_name_fun);
			}
			return ss.str();
		}

	};

	///Prints polynomial
	template<typename container_t>
	std::ostream& operator<<(std::ostream& os, const Polynomial<container_t>& a) {
		os << a.print();
		return os;
	}

}
