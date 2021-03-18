#pragma once
#include <map>
#include <unordered_map>
#include "General.hpp"

///@file
///@brief Contains the class of polynomials in multiple variables.

namespace Symmetric_Polynomials
{

	////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief			The default container for polynomials using the STL map or unordered map
	/// @tparam scl_t	The scalar (coefficient type) (eg \c int,float)
	/// @tparam exp_t	The exponent type eg Standard_Variables, Half_Idempotent_Variables
	/// @tparam ordered Whether to use \c std::map or \c std::unordered_map.
	///////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename scl_t, typename exp_t, bool ordered = 1>
	struct default_container;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief				Polynomials in multiple variables
	///	@tparam container_t Specifies the container for the monomials in the polynomial eg \ref default_container<int, Standard_Variables<>, 1> or \ref default_container<float, Half_Idempotent_Variables<>, 0>
	///	@par				Requirements from container_t
	///						The container must map pairs (degree,exponent) to scalars: \f$(d,e)\mapsto c\f$ represents the monomial \f$cx_1^{e_1}\cdots x_n^{e_n}\f$ of degree \f$d\f$.
	///						Specifically:
	///	- 					\c container_t must have a member \c static \c constexpr \c bool \c ordered that specifies whether it keeps the monomials in increasing order (yes if \c ordered==1, no otherwise)
	///	- 					If \c ordered==1, \c container_t should be equivalent to \c std::map<std::pair<typename exp_t::deg_t, exp_t>, scl_t>
	///	- 					If \c ordered==0, \c container_t should be equivalent to \c std::unordered_map<std::pair<typename exp_t::deg_t, exp_t>, scl_t, hasher>
	///	-					In particular \c container_t must also have:
	///						\c operator[](), \c const_iterator with \c begin(),end() methods, and methods \c insert() and \c erase().
	///	-					Finally, it must have the typedefs \c mapped_type (the scalar), \c key_type::first_type (the degree), \c key_type::second_type (the exponent).
	///	@par				Requirements from the exponent type \c exp_t (i.e. \c key_type::second_type):
	///						The exponent must have functionality similar to \c Standard_Variables or \c Half_Idempotent_Variables
	///						Specifically:
	///	- 					It must have a typedef \c deg_t that represents the degree type (eg \c int,uint64_t)
	///	- 					It must have basic vector functionality (constructor that takes \c int \c n and produces exponent of 0's with that number of variables \c n, operator []...)
	///	- 					It should have an operator + to be used in the product of monomials and an operator - to be used in the division of monomials (if products/divisions need to be used).
	///	- 					If \c ordered==1 it needs to have an operator < (implementing a lexicographic ordering) and otherwise it needs to have an operator () implementing a hash function
	///	- 					It may also optionally have a \c deg_t \c degree() function that compute the degree of the exponent
	///	-					It may also optionally have a \c std::string \c static \c name(int,int) function that prints the names of the variables (first parameter is the index of the variable, second is the total number of variables).
	///	- 					If \c exp_t does not have a \c degree() function then the user needs to provide the dimensions of the variables through a \c deg_t*
	///	- 					If \c exp_t does not have a \c static \c name(int,int) function then the user needs to provide the names of the variables through a \c std::string*
	///	@todo				Use concepts (C++20) to express these requirements
	///	@todo				Improve multiplication \c operator*=()
	/// @note				We save the degree in order to speed up the \c highest_term() computation
	template <typename container_t>
	class Polynomial
	{
	public:
		typedef typename container_t::mapped_type scl_t;		   ///<The scalar type eg \c int
		typedef typename container_t::key_type::second_type exp_t; ///<The exponent type eg \c Standard_Variables
		typedef typename container_t::key_type::first_type deg_t;  ///<The degree type eg \c uint64_t

		///	@brief	Default constructor (empty polynomial)
		Polynomial() = default;

		///	@brief			Constructs empty polynomial with certain number of variables and optionally dimensions/names
		///	@param num_var 	The number of variables
		///	@param dim_var 	Pointer to the dimensions of the variables (used only when \c exp_t does not implement method \c deg_t \c degree())
		///	@param name_var	Pointer to the names of the variables (used only when \c exp_t does not implement \c std::string \c static \c name(int,int))
		/// @attention		The \p dim_var and \p name_var need not be provided (as long as \c exp_t implements the methods above), but if they are,
		///					the user must ensure the pointers are not invalidated
		Polynomial(int num_var, const deg_t *dim_var = nullptr, const std::string *name_var = nullptr);

		///	@brief			Constructs polynomial with a single monomial term
		///	@param exp 		The exponent of the monomial
		///	@param coeff 	The coefficient of the monomial
		///	@param dim_var 	Pointer to the dimensions of the variables (used only when \c exp_t does not implement method \c deg_t \c degree())
		///	@param name_var	Pointer to the names of the variables (used only when \c exp_t does not implement \c std::string \c static \c name(int,int))
		/// @attention		The \p dim_var and \p name_var need not be provided (as long as \c exp_t implements the methods above), but if they are,
		///					the user must ensure the pointers are not invalidated
		Polynomial(const exp_t &exp, scl_t coeff, const deg_t *dim_var = nullptr, const std::string *name_var = nullptr);

		///	@brief			Constructs constant polynomial
		///	@param num_var	The number of variables
		///	@param coeff 	The coefficient of the monomial
		///	@param dim_var 	Pointer to the dimensions of the variables (used only when \c exp_t does not implement method \c deg_t \c degree())
		///	@param name_var	Pointer to the names of the variables (used only when \c exp_t does not implement \c std::string \c static \c name(int,int))
		/// @attention		The \p dim_var and \p name_var need not be provided (as long as \c exp_t implements the methods above), but if they are,
		///					the user must ensure the pointers are not invalidated
		Polynomial(int num_var, scl_t coeff, const deg_t *dim_var = nullptr, const std::string *name_var = nullptr);

		///	@brief	Returns the number of variables of the polynomial
		///	@return The number of variables of the polynomial
		auto number_of_variables() const;

		///	@brief		Returns number of monomials
		///	@attention	Does nothing if the container used does not have a reserve function
		///	@param n	The amount of expected monomials
		void reserve(size_t n);

		///	@brief	Returns the number of monomials of the polynomial
		///	@return The number of of monomials of the polynomial
		auto number_of_monomials() const;

		///	@brief		Clears all monomials of polynomial (equivalent to setting polynomial to 0)
		///	@attention	Does not affect number of variables
		void clear();

		///	@brief				Inserts monomial given exponent and coefficient
		///	@param exponent 	The exponent of the monomial.
		///	@param scalar 		The coefficient of the monomial.
		///	@attention			If a monomial with same exponent already exists in the polynomial then insert does nothing
		///	@warning			It is the user's responsibility to make sure that the coefficient is nonzero
		void insert(const exp_t &exponent, scl_t scalar);

		///	@brief		Constant iterator through the monomials of the polynomial
		///	@warning	The monomials are traversed in increasing order only for ordered containers
		class const_iterator
		{
		public:
			auto coeff() const;					   ///<Returns the coefficient of the monomial
			const auto &exponent() const;		   ///<Returns the exponent of the monomial
			auto degree() const;				   ///<Returns the degree of the monomial
			const_iterator &operator++();		   ///<Increments iterator
			bool operator==(const_iterator) const; ///<Equality of iterators
			bool operator!=(const_iterator) const; ///<Inequality of iterators
		private:
			typename container_t::const_iterator it;
			const_iterator(typename container_t::const_iterator);
			friend class Polynomial; ///<Befriend outer class
		};

		const_iterator begin() const;		 ///<Returns \c const_iterator to the first monomial
		const_iterator end() const;			 ///<Returns \c const_iterator to the end
		const_iterator highest_term() const; ///<Returns \c const_iterator to the highest term monomial

		///	@brief 	Replaces \c *this by \c *this+other.
		///	@note 	Efficient, in place
		/// @return	Reference to \c *this
		Polynomial &operator+=(const Polynomial &other);

		///	@brief 	Replaces \c *this by \c *this-other.
		///	@note 	Efficient, in place
		/// @return	Reference to \c *this
		Polynomial &operator-=(const Polynomial &other);

		///	@brief 		Adds polynomials \c *this and \c other
		/// @return	\c *this+other
		Polynomial operator+(const Polynomial &other) const;

		///	@brief 		Subtracts polynomials \c *this and \c other
		/// @return	\c *this-other
		Polynomial operator-(const Polynomial &other) const;

		///	@brief 		Multiplies polynomials \c *this and \c other
		/// @return	\c *this*other
		Polynomial operator*(const Polynomial &other) const;

		///	@brief 	Replaces \c *this by \c *this*other.
		///	@todo 	Could this be done in place?
		/// @return	Reference to \c *this
		Polynomial &operator*=(const Polynomial &other);

		///	@brief 			Replaces \c *this by \c *this*scalar.
		///	@param	scalar	Scalar we multiply with
		///	@note 			Efficient, in place
		/// @return			Reference to \c *this
		Polynomial &operator*=(scl_t scalar);

		///	@brief 					Raises \c *this to given power \p p
		///	@tparam	any_int_type	Any integer type eg \c int,uint64_t
		/// @warning				Does nothing if \p p<0
		///	@attention				Raises \c static_assert if \c any_int_type is not an integer type
		///	@param 		p			Power we raise to
		///	@todo 					Improve compared to multiplying p many times (eg iterating squaring)
		template <typename any_int_type = int>
		Polynomial operator^(any_int_type p) const;

		///	@brief Prints polynomial
		std::string print() const;

		/// @brief Equality of polynomials
		bool operator==(const Polynomial &) const;

		/// @brief Inequality of polynomials
		bool operator!=(const Polynomial &) const;

	private:
		int _number_of_variables;
		container_t data;
		const deg_t *_dimensions;			//if exp_t does not have the appropriate method
		const std::string *_variable_names; //if exp_t does not have the appropriate method

		auto compute_degree(const exp_t &); //Finds degree of given exponent

		template <typename key_t>
		void insert_add_erase(const key_t &, scl_t);

		template <typename fun>
		void print(scl_t, const exp_t &, std::stringstream &, const fun &) const; //Print monomial using given variable names
		template <typename fun>
		std::string print(const fun &) const; //Print polynomial using given variable names
	};

	/// @brief Prints polynomial to output stream
	template <typename container_t>
	std::ostream &operator<<(std::ostream &os, const Polynomial<container_t> &a);
}
#include "impl/Polynomials.ipp"