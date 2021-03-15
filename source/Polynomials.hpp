#pragma once
#include <map>
#include <unordered_map>
#include "General.hpp"

///@file
///@brief Contains the class of polynomials in multiple variables. 

namespace Symmetric_Polynomials {

	////////////////////////////////////////////////////////////////////////
	///The default container for polynomials using the STL map or unordered_map
    //
	///@tparam scl_t The scalar (coefficient type) (eg int,float)
	///@tparam exp_t The exponent type eg Standard_Variables, Half_Idempotent_Variables
	///@tparam ordered Whether to use std::map or std::unordered_map.
	////////////////////////////////////////////////////////////////////////
	template<typename scl_t, typename exp_t, bool ordered = 1>
	struct default_container;

	///////////////////////////////////////////////////////////////////////////////////
	///Polynomials in multiple variables
	//
	///@tparam container_t Specifies the container for the monomials in the polynomial eg default_container<int, Standard_Variables<>, 1> or default_container<float, Half_Idempotent_Variables<>, 0>
	///
	///
	///In essence the container must map pairs (degree,exponent)=\f$(d,e)\f$ to scalars=\f$c\f$ representing the monomial \f$cx_1^{e_1}\cdots x_n^{e_n}\f$ of degree \f$d\f$.
	///
	///
	///Requirements from container_t:
	///
	///- container_t must have a static constexpr bool variable ordered that specifies whether it keeps the monomials in increasing order (yes if ordered=1, no if ordered=0)
	///
	///- If ordered=1, container_t should be equivalent to std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t> 
	///
	///- If ordered=0, container_t should be equivalent to std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hasher>
	///
	///In particular container_t must also have: operator[](), const_iterator with begin(),end() methods, and methods insert() and erase().
	///
	///Finally, it must have the typedefs mapped_type (the scalar), key_type::first_type (the degree), key_type::second_type (the exponent).
	///
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

		typedef typename container_t::mapped_type scl_t; 			///<The scalar type
		typedef typename container_t::key_type::second_type exp_t;  ///<The exponent type
		typedef typename container_t::key_type::first_type deg_t; 	///<The degree type

		///Default constructor (empty polynomial)
		Polynomial() = default;

		///Constructs empty polynomial
		//
		///@param num_var 	The number of variables
		///@param dim_var 	Pointer to the dimensions of the variables (used only when exp_t does not implement method degree())
		///@param name_var	Pointer to the names of the variables (used only when exp_t does not implement static method name(int,int))
		Polynomial(int num_var, const deg_t* dim_var= nullptr, const std::string* name_var = nullptr);

		///Constructs polynomial with a single monomial term
		//
		///@param exp 		The exponent of the monomial
		///@param coeff 	The coefficient of the monomial
		///@param dim_var 	Pointer to the dimensions of the variables (used only when exp_t does not implement method degree())
		///@param name_var	Pointer to the names of the variables (used only when exp_t does not implement static method name(int,int))
		Polynomial(const exp_t& exp, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///Constructs constant polynomial
		//
		///@param num_var	The number of variables
		///@param coeff 	The coefficient of the monomial
		///@param dim_var 	Pointer to the dimensions of the variables (used only when exp_t does not implement method degree())
		///@param name_var	Pointer to the names of the variables (used only when exp_t does not implement static method name(int,int))
		Polynomial(int num_var, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		auto number_of_variables() const; 	///<Returns the numbers of variables of the polynomial
		void reserve(size_t n);				///<Reserve number of monomials (does nothing if the container used does not have a reserve function)
		auto number_of_monomials() const;	///<Returns number of monomials in polynomial
		void clear();						///<Clears all data of polynomial apart from number of variables

		///Inserts monomial given exponent and coefficient
		//
		///@param exponent 	The exponent of the monomial. If a monomial with same exponent already exists in the polynomial then insert does nothing
		///@param scalar 	The coefficient of the monomial. It is the user's responsibility to make sure that the coefficient is nonzero
		void insert(const exp_t& exponent, scl_t scalar);

		///Constant iterator going through the monomials of the polynomial in increasing order
		class const_iterator {
		public:
			auto coeff() const;						///<Returns the coefficient of the monomial
			const auto& exponent() const;			///<Returns the exponent of the monomial
			auto degree() const;					///<Returns the degree of the monomial
			const_iterator& operator ++();			///<Increments iterator
			bool operator ==(const_iterator) const;	///<Equality of iterators
			bool operator !=(const_iterator) const;	///<Inequality of iterators
		private:
			typename container_t::const_iterator it;
			const_iterator(typename container_t::const_iterator);
			
			///Befriend polynomial
			friend class Polynomial;
		};

		const_iterator begin() const;			///<Returns starting iterator 
		const_iterator end() const;				///<Returns ending iterator
		const_iterator highest_term() const;	///<Returns iterator to the highest term monomial 

		///Adds polynomial on the RHS to polynomial on the LHS (in place)
		Polynomial& operator+=(const Polynomial&);

		///Subtracts polynomial on the RHS from polynomial on the LHS (in place)
		Polynomial& operator-=(const Polynomial&);

		///Adds polynomials
		Polynomial operator+(const Polynomial&) const;

		///Subtracts polynomials
		Polynomial operator-(const Polynomial&) const;

		///Multiplies polynomials
		Polynomial operator* (const Polynomial&) const;

		///Multiplies polynomial on the LHS with polynomial on the RHS (not in place)
		Polynomial& operator*= (const Polynomial&);

		///Multiplies polynomial on the LHS with scalar coefficient on the RHS (in place)
		Polynomial& operator*= (scl_t);

		///Raises polynomial to given power (does nothing if given power is <0)
		template<typename any_int_type = int>
		Polynomial operator^(any_int_type) const;

		///Prints polynomial
		std::string print() const;

		///Equality of polynomials
		bool operator==(const Polynomial&) const;

		///Inequality of polynomials
		bool operator!=(const Polynomial&)  const;

	private:

		int _number_of_variables;
		container_t data;
		const deg_t* _dimensions; //if exp_t does not have the appropriate method
		const std::string* _variable_names; //if exp_t does not have the appropriate method

		auto compute_degree(const exp_t&); 		//Finds degree of given exponent

		template<typename key_t>
		void insert_add_erase(const key_t&, scl_t);

		template<typename fun>
		void print(scl_t, const exp_t&, std::stringstream&, const fun&) const; 	//Print monomial using given variable names
		template<typename fun>
		std::string print(const fun&) const; 		//Print polynomial using given variable names
	};

	///Prints polynomial to output stream
	template<typename container_t>
	std::ostream& operator<<(std::ostream& os, const Polynomial<container_t>& a);
}
#include "impl/Polynomials.ipp"