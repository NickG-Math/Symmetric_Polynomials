#pragma once
#include "General.hpp"
#include "impl/Details.ipp"
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>

///@file
///@brief Contains the class of polynomials in multiple variables.

namespace symmp
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief							Class for polynomials in multiple variables with relations
	/// @tparam _scl					The scalar/coefficient type of the polynomial eg \c float or \c int64_t
	/// @tparam _exp					The variable/exponent type of the polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	/// @tparam _container				The data storage type of the polynomial. 
	///									This should be equivalent to \c std::map if \c container_is_ordered==1 and \c std::unordered_map otherwise
	/// @tparam container_is_ordered	This should be 1 if \c _container is equivalent to \c std::map and 0 if it's equivalent to \c std::unordered_map
	/// @tparam _Args					Any extra optional arguments to pass to \c _container apart from key,value, comparator/hash. Eg an allocator
	/// @par							Requirements from \c exp_t:
	///									The exponent must have functionality similar to \c StandardVariables or \c HalfIdempotentVariables
	///									Specifically:
	///	- 								It must have a typedef \c deg_t that represents the degree type (eg ``` int, uint64_t ```)
	///	- 								It must have basic vector functionality (constructor that takes ```int n``` and produces exponent of 0's with that number of variables \c n, operator []...)
	///	- 								It should have an operator + to be used in the product of monomials (if products need to be used)
	///									and an operator - to be used in the division of monomials (if divisions need to be used).
	///	- 								If \c container_is_ordered==1 it needs to have a ``` bool operator<()  const ```
	///									and otherwise it needs to have a hash function ``` size_t operator()() const ```
	///	- 								Optionally, it may have a ``` deg_t degree() const``` method that computes the degree of the monomial exponent
	///	-								Optionally, it may a ```std::string static name(int,int)``` method that prints the names of the variables 
	///									(first parameter is the index of the variable, second is the total number of variables).
	///	- 								If \c exp_t does not have a \c degree method then the user needs to provide the dimensions of the variables through a \c deg_t* in the constructor
	///	- 								If \c exp_t does not have a \c name method then the user needs to provide the names of the variables through a \c std::string* in the constructor
	///	@todo							Use concepts (C++20) to express these requirements
	template <typename _scl, typename _exp, template<typename...> typename _container, bool container_is_ordered, typename ... _Args>
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
		///	@attention	Does nothing if \c _container does not have a reserve function
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
		///	@note	The monomials are traversed in increasing order only when \c _container is ordered
		class constIterator
		{
		public:
			/// @brief	Returns the scalar coefficient of the monomial the iterator points to
			/// @return The coefficient of the monomial
			auto coeff() const->scl_t;				

			/// @brief	Returns the exponent of the monomial the iterator points to
			/// @return	The exponent of the monomial
			auto exponent() const -> const exp_t&;

			/// @brief	Returns the degree of the monomial the iterator points to
			/// @return	The degree of the monomial
			auto degree() const->deg_t;				

			/// @brief	Increments iterator and points it to the next monomial
			/// @return	Reference to iterator pointing to the next monomial
			auto operator++()->constIterator&;

			/// @brief			Checks if two iterators are not equal (eg comparison to end())
			/// @param	other	The iterator we compare \c *this to
			/// @return			1 if \c *this!=other
			bool operator!=(constIterator other) const;
		private:
			typedef typename implementation_details::polynomial_data_type<_scl, _exp, _container, container_is_ordered, _Args...>::const_iterator it_t;
			it_t it;
			constIterator(it_t);
			friend class Polynomial; ///<Befriend outer class
		};

		/// @brief	\c constIterator to the first monomial
		/// @return \c constIterator pointing to the first monomial
		auto begin() const->constIterator;			

		///	@brief	\c constIterator to the end
		/// @return \c constIterator pointing to just after the last monomial
		/// @note	Only use for comparisons
		auto end() const->constIterator;

		/// @brief	\c constIterator to the highest term monomial
		/// @return \c constIterator pointing to the highest term monomial
		/// @note	If the container is ordered then this just returns the last term (\f$O(1)\f$). Otherwise it finds the highest degree by linear search through the entire polynomial (\f$O(n)\f$)
		auto highest_term() const->constIterator;

		///	@brief 			Addition assignment
		///	@param	other	The polynomial we add to \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		auto operator+=(const Polynomial& other)->Polynomial&;

		///	@brief 			Subtraction assignment
		///	@param	other	The polynomial we subtract from \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		auto operator-=(const Polynomial& other)->Polynomial&;

		///	@brief 			Multiplication assignment
		///	@param	other	The polynomial we multiply with \c *this	
		/// @return			Reference to \c *this
		///	@todo 			Could this be done in place?
		auto operator*=(const Polynomial& other)->Polynomial&;

		///	@brief 			Scalar multiplication assignment
		///	@param	scalar	The scalar we multiply with \c *this	
		/// @return			Reference to \c *this
		///	@note 			Efficient, in place
		auto operator*=(scl_t scalar)->Polynomial&;

		///	@brief 			Addition of polynomials
		///	@param	other	The polynomial we add to \c *this
		/// @return			```(*this)+other```
		auto operator+(const Polynomial& other) const->Polynomial;

		///	@brief 			Subtraction of polynomials
		///	@param	other	The polynomial we subtract from \c *this
		/// @return			```(*this)-other```
		auto operator-(const Polynomial& other) const->Polynomial;

		///	@brief 			Multiplication of polynomials
		///	@param	other	The polynomial we multiply with \c *this
		/// @return			```(*this)*other```
		auto operator*(const Polynomial& other) const->Polynomial;

		///	@brief 			Raises polynomial to integer power
		///	@tparam		T	Any integer type eg \c int,uint64_t
		///	@param 		p	Power we raise \c *this to
		/// @return			```(*this)^p```
		/// @warning		Does nothing if \p p<0
		///	@attention		Raises \c static_assert if \c T is not an integer type
		///	@todo 			Improve implementation (current is multiplication p-many times; maybe square p/2 times instead?)
		template <typename T = int>
		auto operator^(T p) const->Polynomial;

		/// @brief			Checks if two polynomials are equal
		/// @param	other	The polynomial we compare \c *this to
		/// @return			1 if ```*this==other```
		bool operator==(const Polynomial& other) const;

		/// @brief			Checks if two polynomials are not equal
		/// @param	other	The polynomial we compare \c *this to
		/// @return			1 if ```*this!=other```
		bool operator!=(const Polynomial& other) const;

	private:
		typedef typename implementation_details::polynomial_data_type<_scl, _exp, _container, container_is_ordered, _Args...> data_t;
		data_t data;
		const deg_t* _dimensions;			//if exp_t does not have the appropriate method
		const std::string* _variable_names; //if exp_t does not have the appropriate method
		auto compute_degree(const exp_t&)->deg_t; //Finds degree of given exponent
		template <typename key_t>
		void insert_add_erase(const key_t&, scl_t);
		template <typename fun>
		void print(scl_t, const exp_t&, std::ostream&, const fun&) const; //Print monomial using given variable names
		template <typename fun>
		void print(std::ostream& os, const fun&) const; //Print polynomial using given variable names

		template<typename t1, typename t2, template<typename...> typename t3, bool t4, typename ... t5>
		friend std::ostream& operator<<(std::ostream& os, const Polynomial<t1,t2,t3,t4,t5...>& a); ///<Befriending the print to ostream function
	};

	/// @brief							Prints polynomial to output stream
	/// @tparam scl_t					The scalar/coefficient type of the polynomial
	/// @tparam exp_t					The variable/exponent type of the polynomial
	/// @tparam container_t				The data storage type of the polynomial. 
	/// @tparam container_is_ordered	Whether the container is ordered or not
	/// @tparam Args					Any extra optional arguments to pass to \c _container apart from key,value, comparator/hash
	/// @param	os						The output stream
	/// @param	a						The polynomial to be printed
	/// @return							The output stream where the polynomial has been printed
	template <typename scl_t, typename exp_t, template<typename...> typename container_t, bool container_is_ordered, typename ... Args>
	std::ostream& operator<<(std::ostream& os, const Polynomial<scl_t, exp_t, container_t, container_is_ordered, Args...>& a);

	/// @brief			A polynomial whose monomials are stored in increasing order
	/// @tparam _scl	The scalar/coefficient type of the polynomial
	/// @tparam _exp	The variable/exponent type of the Polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	template <typename _scl, typename _exp>
	using OrderedPolynomial = Polynomial<_scl, _exp, std::map, 1>;

	/// @brief			A polynomial whose monomials are not stored in any particular order
	/// @tparam _scl	The scalar/coefficient type of the polynomial
	/// @tparam _exp	The variable/exponent type of the Polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	template <typename _scl, typename _exp>
	using UnorderedPolynomial = Polynomial<_scl, _exp, std::unordered_map, 0>;
}
#include "impl/Polynomials.ipp"