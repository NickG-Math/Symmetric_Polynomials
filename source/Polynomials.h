#pragma once
#include <functional>
#include <map>
#include "General.h"

///@file
///@brief Contains the class of polynomials in multiple variables. 

namespace Symmetric_Polynomials {

	///Class of Polynomials in multiple variables
	//
	///scalar_t is the coefficient type of the polynomial eg int or Rational
	///
	///exponent_t is the type of exponents used in the polynomial. Eg standard<T> or half_idempotent<T>. exponent_t should have a "int degree()" function and a "static std::string name(int)" function.
	///
	///name_override should be turned to 1 when the names of the variables are provided and the "name" function of exponent_t will not be used. By default 0.
	///
	///dimension_override should be turned to 1 when the dimensions of the variables are provided and the "degree" function of exponent_t will not be used. By default 0.
	template<typename scalar_t, typename exponent_t, bool name_override = 0, bool dimension_override = 0>
	class Polynomial {
		std::map<std::pair<int, exponent_t>, scalar_t> data;
		std::vector<std::string> variable_names; //if exponent_t does not have the appropriate method
		std::vector<int> dimensions; //if exponent_t does not have the appropriate method
		int number_of_variables;
	public:
		///Default constructor
		Polynomial() {}

		///Constructor usuing number of variables
		Polynomial(int number_of_variables) : number_of_variables(number_of_variables) {}

		///Constructor using variable names and dimensions (name_override=1 and dimension_override=1)
		Polynomial(const std::vector<std::string>& names, const std::vector<int>& dims) : Polynomial(names.size()) { set_variable_names(names); set_dimensions(dims); }

		///Constructor using variable names (name_override=1)
		Polynomial(const std::vector<std::string>& names) : Polynomial(names.size()) { set_variable_names(names); }

		///Constructor using variable dimensions (dimension_override=1)
		Polynomial(const std::vector<int>& dims) : Polynomial(dims.size()) { set_dimensions(dims); }

		///Constructs monomial given exponent and scalar coefficient
		Polynomial(const exponent_t& exponent, const scalar_t& scalar) : Polynomial(exponent.size()){insert(exponent, scalar);}

		///Constructs monomial given exponent and scalar coefficient and using variable names and dimensions (name_override=1 and dimension_override=1)
		Polynomial(const exponent_t& exponent, const scalar_t& scalar, const std::vector<std::string>& names, const std::vector<int>& dims) : Polynomial(names, dims){insert(exponent, scalar);}

		///Constructs monomial given exponent and scalar coefficient and using variable names (name_override=1)
		Polynomial(const exponent_t& exponent, const scalar_t& scalar, const std::vector<std::string>& names) : Polynomial(names){insert(exponent, scalar);}

		///Constructs monomial given exponent and scalar coefficient and using using variable dimensions (dimension_override=1)
		Polynomial(const exponent_t& exponent, const scalar_t& scalar, const std::vector<int>& dims) : Polynomial(dims) {insert(exponent, scalar);}

		///Returns the numbers of variables of the polynomial
		int get_number_of_variables() const{return number_of_variables;}

		///Returns degree of polynomial
		int degree() const {
			auto last = std::prev(data.end());
			return last->first.first;
		}

		///Clears all data of polynomial apart from number of variables
		void clear() {data.clear();}

		///Inserts monomial given exponent and coefficient
		void insert(const exponent_t& exponent, scalar_t scalar) {
			if constexpr (dimension_override)
				data[std::pair<int, exponent_t>(general_compute_degree(exponent, dimensions), exponent)] = scalar;
			else
				data[std::pair<int, exponent_t>(exponent.degree(), exponent)] = scalar;
		}


		///Constant iterator going through the monomials of the polynomial in increasing order
		class const_iterator {
		public:
			///Returns the coefficient of the monomial
			const scalar_t& coeff() const {
				return it->second;
			}
			///Returns the exponent of the monomial
			const exponent_t& exponent() const {
				return (it->first).second;
			}
			///Returns the degree of the monomial
			int degree() const {
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
			const_iterator(typename std::map<std::pair<int, exponent_t>, scalar_t>::const_iterator it) : it(it) {}
			typename std::map<std::pair<int, exponent_t>, scalar_t>::const_iterator it;
			template<typename, typename, bool, bool>
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
			return std::prev(data.end());
		}

		///Constructs constant polynomial given coefficient and number of variables
		static Polynomial<scalar_t, exponent_t> constant(int coeff, int variables) {
			Polynomial<scalar_t, exponent_t> c(variables);
			exponent_t v = exponent_t(variables);
			c.data[std::pair<int, exponent_t>(0, v)] = coeff;
			return c;
		}

		///Adds polynomial on the RHS to polynomial on the LHS (in place)
		Polynomial<scalar_t, exponent_t>& operator+=(const Polynomial<scalar_t, exponent_t>& b) {
			for (const auto& pair : b.data) {
				data[pair.first] += pair.second;
				if (data[pair.first] == 0)
					data.erase(pair.first);
			}
			return *this;
		}

		///Subtracts polynomial on the RHS from polynomial on the LHS (in place)
		Polynomial<scalar_t, exponent_t>& operator-=(const Polynomial<scalar_t, exponent_t>& b) {
			for (const auto& pair : b.data) {
				data[pair.first] -= pair.second;
				if (data[pair.first] == 0)
					data.erase(pair.first);
			}
			return *this;
		}

		///Adds polynomials
		Polynomial<scalar_t, exponent_t> operator+(const Polynomial<scalar_t, exponent_t>& b) const {
			Polynomial<scalar_t, exponent_t> sum = *this;
			sum += b;
			return sum;
		}

		///Subtracts polynomials
		Polynomial<scalar_t, exponent_t> operator-(const Polynomial<scalar_t, exponent_t>& b) const {
			Polynomial<scalar_t, exponent_t> diff = *this;
			diff -= b;
			return diff;
		}

		///Multiplies polynomials
		Polynomial<scalar_t, exponent_t> operator* (const Polynomial<scalar_t, exponent_t>& b) const {
			Polynomial<scalar_t, exponent_t> product(number_of_variables);
			for (const auto& paira : data) {
				for (const auto& pairb : b.data) {
					auto vsum = std::pair(paira.first.first + pairb.first.first, paira.first.second + pairb.first.second);
					product.data[vsum] += paira.second * pairb.second;
					if (product.data[vsum] == 0)
						product.data.erase(vsum);
				}
			}
			return product;
		}

		///Multiplies polynomial on the LHS with polynomial on the RHS (not in place)
		Polynomial<scalar_t, exponent_t>& operator*= (const Polynomial<scalar_t, exponent_t>& b) {
			*this = operator*(b);
			return *this;
		}

		///Multiplies polynomial on the LHS with scalar coefficient on the RHS (in place)
		Polynomial<scalar_t, exponent_t>& operator*= (const scalar_t& coeff) {
			if (coeff == 0)
				clear();
			else if (coeff!=1){
				for (const auto& paira : data)
					data[paira.first] *= coeff;
			}
			return *this;
		}

		///Raises polynomial to given power
		Polynomial<scalar_t, exponent_t> operator^(int p) const {
			if (p == 0)
				return constant(1, number_of_variables);
			if (p == 1)
				return *this;
			auto prod = *this;
			for (int k = 1; k < p; k++)
				prod *= *this;
			return prod;
		}

		///Equality of polynomials
		bool operator==(const Polynomial<scalar_t, exponent_t>& b) const {
			return (data == b.data);
		}

		///Inequality of polynomials
		bool operator!=(const Polynomial<scalar_t, exponent_t>& b) const {
			return !(*this == b);
		}

		///Print polynomial
		std::string print() const {
			if constexpr (name_override)
				return print([=](int i) {return variable_names[i];});
			else
				return print(exponent_t::name);
		}

	private:

		void set_variable_names(const std::vector<std::string>& names) {
			if constexpr (!name_override)
				static_assert("The template parameter name_override must be set to 1 to allow custom names");
			variable_names = names;
		}

		void set_dimensions(const std::vector<int>& dims) {
			if constexpr (!dimension_override)
				static_assert("The template parameter dimension_override must be set to 1 to allow custom dimensions");
			dimensions = dims;
		}


		///Print monomial using given variable names
		void print(const scalar_t& coeff, const exponent_t& exponent, std::stringstream& ss, const std::function<std::string(int)>& variable_names) const {
			bool havestar = 0;
			if (coeff != 1) {
				ss << coeff;
				havestar = 1;
			}
			bool completelyzero = 1;
			for (decltype(exponent.size()) i = 0; i < exponent.size(); i++) {
				if (exponent[i] != 0) {
					completelyzero = 0;
					if (havestar)
						ss << "*";
					havestar = 1;
					if (exponent[i] > 1) {
						ss << variable_names(i) << "^" << (int)exponent[i];
					}
					else if (exponent[i] == 1)
						ss << variable_names(i);
				}
			}
			if (completelyzero && coeff == 1)
				ss << "1";
		}


		///Print polynomial using given variable names
		std::string print(const std::function<std::string(int)>& variable_name_fun) const {
			std::stringstream ss;
			auto it = data.begin();
			print(it->second, it->first.second, ss, variable_name_fun);
			for (it++; it != data.end(); it++) {
				ss << " + ";
				print(it->second, it->first.second, ss, variable_name_fun);
			}
			return ss.str();
		}
	
	};

	///Prints polynomial
	template<typename scalar_t, typename exponent_t, bool t, bool w>
	std::ostream& operator<<(std::ostream& os, const Polynomial<scalar_t, exponent_t, t, w>& a) {
		os << a.print();
		return os;
	}

	///The standard variables \f$x_i\f$ in a polynomial, with \f$|x_i|=1\f$ and no relations.
	template<typename T>
	struct Standard_Variables : public std::vector<T> {
		using std::vector<T>::vector;
		///Degree of exponent
		int degree() const {
			return sum(*this);
		}
		///Name of variable
		static std::string name(int i) {
			return "x_" + std::to_string(i + 1);
		}
		///Addition of exponents
		Standard_Variables<T> operator+ (const Standard_Variables<T>& other) const {
			Standard_Variables<T> v(*this);
			for (int i = 0; i < this->size(); i++)
				v[i] += other[i];
			return v;
		}
	};


}
