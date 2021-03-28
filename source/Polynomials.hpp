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
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief			The base of all monomial containers
	/// @tparam	 exp_t	The variable/exponent type of the monomials
	/// @attention		\p exp_t must have functionality similar to \c StandardVariables, \c ElementarySymmetricVariables,
	///					\c HalfIdempotentVariables or \c TwistedChernVariables. Specifically:
	///	- 				It must have a typedef \c deg_t that represents the degree type (eg ``` int, size_t ```)
	///	- 				It must have basic vector functionality (easiest way is to inherit from ```std::vector<size_t>```)
	///	- 				It should have operators + and - to be used in the product/division of monomials
	/// -				Optionally, it may have a ``` deg_t degree() const``` method that computes the degree of the monomial exponent
	///	@warning		If \c exp_t does not have a \c degree method then the dimensions of the variables must be provided through 
	///					a valid ``` deg_t* dimensions``` in the constructor
	///	@todo			Use concepts (C++20) to express these requirements
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class exp_t>
	class BaseContainer {
		typedef typename exp_t::deg_t deg_t;
	protected:

		/// @brief				Constructor given pointer to dimensions of the variables
		/// @param dimensions	Used only when \c exp_t does not implement method ```deg_t degree() const```
		BaseContainer(const deg_t* dimensions = nullptr);

		/// @brief		Computes degree of given exponent
		/// @note		Calls ```typename exp_t::degree()``` if that method exists, otherwise uses \ref dimensions
		///	@warning	If neither \ref dimensions is set nor ```typename exp_t::degree()``` exists then this function calls \c abort()
		deg_t compute_degree(const exp_t& exponent) const;

		const deg_t* dimensions;	///<The pointer to the dimensions of the variables
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief				The default ordered/unordered monomial container
	/// @details			Monomials are stored as key-value pairs where each key consists of degree+exponent and 
	///						the value is the coefficient of the monomial.
	///	@note				We store the degree for improved performance.
	/// @tparam _scl		The scalar/coefficient type of the monomials eg \c float or \c int64_t
	/// @tparam _exp		The variable/exponent type of the monomials eg \c StandardVariables or \c HalfIdempotentVariables
	/// @tparam _cnt		The underlying container: this should be equivalent to \c std::map if \c _ord==1 and \c std::unordered_map otherwise
	/// @tparam _ord		This should be 1 if \c _cnt is equivalent to \c std::map and 0 if it's equivalent to \c std::unordered_map
	/// @tparam _arg		Any extra optional arguments to pass to \c _cnt apart from key,value, comparator/hash 
	///						(the latter two are provided in \c _exp) eg an allocator
	///	@attention			In addition to the requirements from BaseContainer, 
	///						if \c _ord==1 then \c exp_t needs to have a comparator ``` bool operator<(const _exp&)  const ```
	///						and otherwise \c exp_t needs to have a hash function ``` size_t operator()() const ```
	///	@todo				Use concepts (C++20) to express these requirements
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
	class DefaultContainer : public BaseContainer<_exp>, implementation_details::container_wrapper<_scl, _exp, _cnt, _ord, _arg...> {
	private:
		typedef typename implementation_details::container_wrapper<_scl, _exp, _cnt, _ord, _arg...> data_t;
	protected:
		typedef _scl scl_t;					///<The scalar (coefficient) type eg \c int64_t
		typedef _exp exp_t;					///<The exponent (variable) type eg \c StandardVariables or \c HalfIdempotentVariables
		typedef typename _exp::deg_t deg_t;	///<The degree type eg \c size_t
	public:
		/// @brief	Returns number of monomials in the polynomial
		size_t number_of_monomials() const;

		/// @brief Equality of polynomials
		bool operator==(const DefaultContainer&) const;

		/// @brief Inequality of polynomials
		bool operator!=(const DefaultContainer& b) const;

		/// @brief			Reserve number of monomials in polynomial
		///	@param		n	Number to reserve
		/// @attention		Does nothing if ```_ord==1```
		void reserve(size_t n);

		/// @brief			Inserts monomial in polynomial.
		///	@param	exp		The exponent of the monomial.
		///	@param	coeff	The coefficient of the monomial.
		/// @warning		It is the user's responsibility to make sure that all exponents have the same length(number of variables),
		///					the exponent being inserted is not already there and that \c coeff is not 0
		void insert(const exp_t& exp, scl_t coeff);

		/// @brief Constant iterator traversing the monomials of a polynomial
		class ConstIterator : public data_t::const_iterator {
		public:
			scl_t coeff() const;			///<The coefficient of the monomial
			const exp_t& exponent() const;	///<The exponent of the monomial
			deg_t degree() const;			///<The degree of the monomial
		private:
			ConstIterator(typename data_t::const_iterator);
			friend class DefaultContainer; ///<Befriending outer class
		};
		ConstIterator begin() const;	///<ConstIterator to the first monomial
		ConstIterator end() const;		///<ConstIterator to just after the final monomial

		/// @brief	ConstIterator to the highest term monomial
		///	@note	If the container is ordered then this just returns the last term in \f$O(1)\f$ time.\n
		///			Otherwise it finds the highest degree by linear search through the entire polynomial in \f$O(n)\f$ time.
		ConstIterator highest_term() const;

	protected:
		using BaseContainer<exp_t>::BaseContainer;
		/// @brief Non const iterator traversing the monomials of a polynomial
		class Iterator : public data_t::iterator {
		public:
			scl_t& coeff();	///<Reference to the coefficient of the monomial
		private:
			Iterator(typename data_t::iterator);
			friend class DefaultContainer; ///<Befriending outer class
		};

		Iterator begin();	///<Iterator to the first monomial
		Iterator end();		///<Iterator to just after the final monomial

		/// @brief		Adds given monomial to polynomial
		/// @param kvp	The key-value pair representing the monomial
		void add(const std::pair<const std::pair<deg_t, exp_t>, scl_t>& kvp);

		/// @brief		Subtracts given monomial from polynomial
		/// @param kvp	The key-value pair representing the monomial
		void subtract(const std::pair<const std::pair<deg_t, exp_t>, scl_t>& kvp);

		/// @brief		Multiplies the two given monomials and then adds the product to polynomial
		/// @param kvp1	The key-value pair representing the first monomial
		/// @param kvp2	The key-value pair representing the second monomial
		void multiply_add(const std::pair<const std::pair<deg_t, exp_t>, scl_t >& kvp1, const std::pair<const std::pair <deg_t, exp_t>, scl_t >& kvp2);
	private:
		void add(const std::pair<deg_t, exp_t> key, scl_t value);
	};


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief				Class for polynomials in multiple variables with relations
	/// @tparam container_t	The data storage type/ monomial container to use eg \c DefaultContainer . \n
	///						The requirements from \c container_t are all non-private methods and typedefs in \c DefaultContainer .
	///	@attention			If \c exp_t does not have a ```std::string static name(int,int)``` method,
	///						then in order to print the polynomial, the names of the variables must be manually provided 
	///						through a valid \c std::string* in the constructor
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename container_t>
	class Polynomial : public container_t
	{
	public:
		typedef typename container_t::scl_t scl_t;	///<The scalar (coefficient) type eg \c int64_t
		typedef typename container_t::exp_t exp_t;	///<The exponent (variable) type eg \c StandardVariables or \c HalfIdempotentVariables
		typedef typename container_t::deg_t deg_t;	///<The degree type eg \c size_t

		///	@brief				Constructs zero polynomial
		///	@param	dim_var		Pointer to the dimensions of the variables; used when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param	name_var	Pointer to the names of the variables; used when \c exp_t does not implement ``` std::string  static name(int,int)```
		Polynomial(const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief				Constructs polynomial with a single nonzero monomial term
		///	@param	exp 		The exponent of the monomial
		///	@param	coeff 		The coefficient of the monomial
		///	@param	dim_var 	Pointer to the dimensions of the variables; used when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param  name_var	Pointer to the names of the variables; used when \c exp_t does not implement ``` std::string  static name(int,int)```
		/// @warning			It is the user's responsibility to make sure \c coeff!=0
		Polynomial(const exp_t& exp, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief				Constructs constant nonzero polynomial
		///	@param	num_var		The number of variables
		///	@param	coeff 		The coefficient of the monomial
		///	@param	dim_var 	Pointer to the dimensions of the variables; used when \c exp_t does not implement method ``` deg_t degree() const ```
		///	@param  name_var	Pointer to the names of the variables; used when \c exp_t does not implement ``` std::string  static name(int,int)```
		/// @warning			It is the user's responsibility to make sure \c coeff!=0
		Polynomial(int num_var, scl_t coeff, const deg_t* dim_var = nullptr, const std::string* name_var = nullptr);

		///	@brief		Returns the number of variables of the polynomial
		///	@return		The number of variables of \c *this
		///	@warning	May only be used on nonempty polynomials
		size_t number_of_variables() const;

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

	private:
		const std::string* variable_names; //if exp_t does not have the appropriate method
		template <class fun>
		void print(scl_t, const exp_t&, std::ostream&, const fun&) const; //Print monomial using given variable names
		template <class fun>
		void print(std::ostream& os, const fun&) const; //Print polynomial using given variable names
		template<class cont>
		friend std::ostream& operator<< (std::ostream&, const Polynomial<cont>&); ///<Befriending the print to ostream function
	};

	/// @brief				Prints polynomial to output stream
	/// @tparam container_t	The data store type (container) of the polynomial eg \c DefaultContainer
	/// @param	os			The output stream
	/// @param	a			The polynomial to be printed
	/// @return				The output stream where the polynomial has been printed
	template <class container_t>
	std::ostream& operator<<(std::ostream& os, const Polynomial<container_t>& a);

	/// @brief			Polynomial using the default containers \c std::map or \c std::unordered_map
	/// @tparam _scl	The scalar/coefficient type of the polynomial
	/// @tparam _exp	The variable/exponent type of the Polynomial eg \c StandardVariables or \c HalfIdempotentVariables
	/// @tparam _ord	If ```_ord==1``` then \c std::map is used as the container for the polynomial; otherwise \c std::unordered_map is used
	template <class _scl, class _exp, bool _ord = 1>
	using Poly = Polynomial<std::conditional_t<_ord, DefaultContainer<_scl, _exp, std::map, 1>, DefaultContainer<_scl, _exp, std::unordered_map, 0>>>;

}
#include "impl/Polynomials.ipp"