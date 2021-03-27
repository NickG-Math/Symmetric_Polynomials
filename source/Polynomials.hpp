#pragma once
#include "General.hpp"
#include "impl/Details.ipp"
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>

///	@file
///	@brief Contains the class of polynomials in multiple variables.

namespace symmp
{

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief				Class for polynomials in multiple variables with relations
	/// @tparam _scl		The scalar/coefficient type of the polynomial eg \c float or \c int64_t
	/// @tparam _exp		The variable/exponent type of the polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	/// @tparam _cnt		The data storage type (container) of the polynomial. 
	///						This should be equivalent to \c std::map if \c _ord==1 and \c std::unordered_map otherwise
	/// @tparam _ord		This should be 1 if \c _cnt is equivalent to \c std::map and 0 if it's equivalent to \c std::unordered_map
	/// @tparam _arg		Any extra optional arguments to pass to \c _cnt apart from key,value, comparator/hash (the latter two are provided in \c _exp). For example an allocator
	/// @par				Requirements from \c exp_t:
	///						The exponent must have functionality similar to \c StandardVariables, \c ElementarySymmetricVariables, \c HalfIdempotentVariables or \c TwistedChernVariables
	///						Specifically:
	///	- 					It must have a typedef \c deg_t that represents the degree type (eg ``` int, uint64_t ```)
	///	- 					It must have basic vector functionality (constructor that takes ```int n``` and produces exponent of 0's with that number of variables \c n, ```operator []```...)
	///	- 					It should have an operator + to be used in the product of monomials (if products need to be used)
	///						and an operator - to be used in the division of monomials (if divisions need to be used).
	///	- 					If \c _ord==1 it needs to have a comparator ``` bool operator<(const _exp&)  const ```
	///						and otherwise it needs to have a hash function ``` size_t operator()() const ```
	///	- 					Optionally, it may have a ``` deg_t degree() const``` method that computes the degree of the monomial exponent
	///	-					Optionally, it may a ```std::string static name(int var,int num_var)``` method that prints the name of each variable
	///	- 					If \c exp_t does not have a \c degree method then the user needs to provide the dimensions of the variables through a \c deg_t* in the constructor
	///	- 					If \c exp_t does not have a \c name method then the user needs to provide the names of the variables through a \c std::string* in the constructor
	///	@todo				Use concepts (C++20) to express these requirements
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	class Polynomial
	{
	public:
		///	@brief	The scalar/coefficient type eg \c int
		typedef _scl scl_t;					
		///	@brief	The variable/exponent type eg \c StandardVariables
		typedef _exp exp_t;				
		///	@brief	The degree type eg \c uint64_t
		typedef typename _exp::deg_t deg_t;

		///	@brief				Constructs zero polynomial
		///	@param	dim_var		Pointer to the dimensions of the variables; used only when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param	name_var	Pointer to the names of the variables; used only when \c exp_t does not implement ``` std::string  static name(int,int)```
		Polynomial(const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief				Constructs polynomial with a single nonzero monomial term
		///	@param	exp 		The exponent of the monomial
		///	@param	coeff 		The coefficient of the monomial
		///	@param	dim_var 	Pointer to the dimensions of the variables; used only when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param  name_var	Pointer to the names of the variables; used only when \c exp_t does not implement ``` std::string  static name(int,int)```
		/// @warning			It is the user's responsibility to make sure \c coeff!=0
		Polynomial(const exp_t& exp, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief				Constructs constant nonzero polynomial
		///	@param	num_var		The number of variables
		///	@param	coeff 		The coefficient of the monomial
		///	@param	dim_var 	Pointer to the dimensions of the variables; used only when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param  name_var	Pointer to the names of the variables; used only when \c exp_t does not implement ``` std::string  static name(int,int)```
		/// @warning			It is the user's responsibility to make sure \c coeff!=0
		Polynomial(int num_var, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief		Reserves number of monomials
		///	@param	n	The amount of expected monomials
		///	@attention	Has no effect if \c _cnt does not have a reserve function
		void reserve(size_t n);

		///	@brief		Returns the number of variables of the polynomial
		///	@return		The number of variables of \c *this
		///	@warning	May only be used on nonempty polynomials
		size_t number_of_variables() const;

		///	@brief	Returns the number of monomials of the polynomial
		///	@return The number of monomials of \c *this
		size_t number_of_monomials() const;

		///	@brief				Inserts monomial in polynomial
		///	@param	exponent 	The exponent of the monomial.
		///	@param	coeff 		The coefficient of the monomial.
		///	@attention			If a monomial with same exponent already exists in the polynomial then insert does nothing
		///	@warning			It is the user's responsibility to make sure that the coefficient is nonzero and that all exponents have the same size (number of variables)
		void insert(const exp_t& exponent, scl_t coeff);

		///	@brief	Constant iterator through the monomials of the polynomial
		///	@note	The monomials are traversed in increasing order only when \c _cnt is ordered
		class constIterator
		{
		public:
			/// @brief	Returns the scalar coefficient of the monomial the iterator points to
			/// @return The coefficient of the monomial
			scl_t coeff() const;				

			/// @brief	Returns the exponent of the monomial the iterator points to
			/// @return	The exponent of the monomial
			const exp_t& exponent() const;

			/// @brief	Returns the degree of the monomial the iterator points to
			/// @return	The degree of the monomial
			deg_t degree() const;				

			/// @brief	Increments iterator and points it to the next monomial
			/// @return	Reference to iterator pointing to the next monomial
			constIterator& operator++();

			/// @brief			Checks if two iterators are not equal (eg comparison to end())
			/// @param	other	The iterator we compare \c *this to
			/// @return			1 if \c *this!=other
			bool operator!=(constIterator other) const;
		private:
			typedef typename implementation_details::polynomial_data_type<_scl, _exp, _cnt, _ord, _arg...>::const_iterator it_t;
			it_t it;
			constIterator(it_t);
			friend class Polynomial; ///<Befriend outer class
		};

		/// @brief	\c constIterator to the beginning
		/// @return \c constIterator pointing to the first monomial
		constIterator begin() const;			

		///	@brief	\c constIterator to the end
		/// @return \c constIterator pointing to just after the last monomial
		/// @note	Only use with ```!=```
		constIterator end() const;

		/// @brief	\c constIterator to the highest term monomial
		/// @return \c constIterator pointing to the highest term monomial
		/// @note	If the container is ordered then this just returns the last term (\f$O(1)\f$). Otherwise it finds the highest degree by linear search through the entire polynomial (\f$O(n)\f$)
		constIterator highest_term() const;

		///	@brief 			Addition assignment
		///	@param	other	The polynomial we add to \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		Polynomial& operator+=(const Polynomial& other);

		///	@brief 			Subtraction assignment
		///	@param	other	The polynomial we subtract from \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		Polynomial& operator-=(const Polynomial& other);

		///	@brief 			Multiplication assignment
		///	@param	other	The polynomial we multiply with \c *this	
		/// @return			Reference to \c *this
		///	@todo 			Could this be done in place?
		Polynomial& operator*=(const Polynomial& other);

		///	@brief 			Scalar multiplication assignment
		///	@param	scalar	The scalar we multiply with \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		Polynomial& operator*=(scl_t scalar);

		///	@brief 			Addition of polynomials
		///	@param	other	The polynomial we add to \c *this
		/// @return			```(*this)+other```
		Polynomial operator+(const Polynomial& other) const;

		///	@brief 			Subtraction of polynomials
		///	@param	other	The polynomial we subtract from \c *this
		/// @return			```(*this)-other```
		Polynomial operator-(const Polynomial& other) const;

		///	@brief 			Multiplication of polynomials
		///	@param	other	The polynomial we multiply with \c *this
		/// @return			```(*this)*other```
		Polynomial operator*(const Polynomial& other) const;

		///	@brief 			Raises polynomial to integer power
		///	@tparam		T	Any integer type eg ``` int, uint64_t```
		///	@param 		p	Power we raise \c *this to
		/// @return			```(*this)^p```
		/// @warning		Does nothing if \p p<0
		///	@attention		Raises \c static_assert if \c T is not an integer type
		///	@todo 			Improve implementation (current is multiplication p-many times; maybe square p/2 times instead?)
		template <class T = int>
		Polynomial operator^(T p) const;

		/// @brief			Checks if two polynomials are equal
		/// @param	other	The polynomial we compare \c *this to
		/// @return			1 if ```*this==other```
		bool operator==(const Polynomial& other) const;

		/// @brief			Checks if two polynomials are not equal
		/// @param	other	The polynomial we compare \c *this to
		/// @return			1 if ```*this!=other```
		bool operator!=(const Polynomial& other) const;

	private:
		typedef typename implementation_details::polynomial_data_type<_scl, _exp, _cnt, _ord, _arg...> data_t;
		data_t data;
		const deg_t* _dimensions;			//if exp_t does not have the appropriate method
		const std::string* _variable_names; //if exp_t does not have the appropriate method
		auto compute_degree(const exp_t&)->deg_t; //Finds degree of given exponent
		template <class key_t>
		void insert_add_erase(const key_t&, scl_t);
		template <class fun>
		void print(scl_t, const exp_t&, std::ostream&, const fun&) const; //Print monomial using given variable names
		template <class fun>
		void print(std::ostream& os, const fun&) const; //Print polynomial using given variable names

		template<class t1, class t2, template<class...> class t3, bool t4, class ... t5>
		friend std::ostream& operator<< (std::ostream&, const Polynomial<t1, t2, t3, t4, t5...>&); ///<Befriending the print to ostream function
	};

	/// @brief			Prints polynomial to output stream
	/// @tparam _scl	The scalar/coefficient type of the polynomial
	/// @tparam _exp	The variable/exponent type of the polynomial
	/// @tparam _cnt	The data storage type (container) of the polynomial. 
	/// @tparam _ord	Whether the container is ordered or not
	/// @tparam _arg	Any extra optional arguments to pass to the container apart from key,value, comparator/hash
	/// @param	os		The output stream
	/// @param	a		The polynomial to be printed
	/// @return			The output stream where the polynomial has been printed
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	std::ostream& operator<<(std::ostream& os, const Polynomial<_scl, _exp, _cnt, _ord, _arg...>& a);

	/// @brief			Polynomial using the default containers \c std::map or \c std::unordered_map
	/// @tparam _scl	The scalar/coefficient type of the polynomial
	/// @tparam _exp	The variable/exponent type of the Polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	/// @tparam _ord	If ```_ord==1``` then \c std::map is used as the container for the polynomial; otherwise \c std::unordered_map is used
	template <class _scl, class _exp, bool _ord=1>
	using Poly = std::conditional_t<_ord, Polynomial<_scl, _exp, std::map, 1>, Polynomial<_scl, _exp, std::unordered_map, 0>>;

}
#include "impl/Polynomials.ipp"