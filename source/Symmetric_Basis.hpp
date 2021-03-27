#pragma once
#include "Polynomials.hpp"
#include "Generators.hpp"

/////////////////////////////////////////////////////////////////////////
///	@file
///	@brief 		Contains the methods and classes for generatic symmetric polynomials
///	@details	The goal is to write any symmetric polynomial with no relations in terms of elementary symmetric polynomials
///				The general interface for doing this can be generalized to subrings of polynomial rings with relations
/////////////////////////////////////////////////////////////////////////

namespace symmp
{

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief 				The standard variables \f$x_i\f$ in a polynomial, with \f$|x_i|=1\f$ and no relations.
	///	@details			A monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ is stored as the vector \f$[a_1,...,a_n]\f$
	///	@tparam		T 		The (integral) value type of the exponent vector.
	///	@tparam		_deg	The (integral) value type used in the degree function.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class T = int64_t, class _deg = int64_t>
	struct StandardVariables : public std::vector<T>
	{
		using std::vector<T>::vector;
		typedef _deg deg_t; ///<Degree typedef

		///	@brief	Computes degree of monomial on standard variables \f$x_i\f$
		///	@return Degree \f$\sum_ia_i\f$ for monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ (```*this```=\f$[a_1,...,a_n]\f$)
		deg_t degree() const;

		///	@brief		Returns the names of the standard variables \f$x_i\f$
		///	@return  	\c "x_i"
		///	@param	i 	The variable index
		///	@param 	n 	The number variables
		static std::string name(int i, int n);

		///	@brief		Multiplies monomials by adding their exponents
		///	@param 	b 	The exponent \f$[b_1,...,b_n]\f$ we add to ```*this```=\f$[a_1,...,a_n]\f$
		///	@return 	Degree \f$[a_1+b_1,...,a_n+b_n]\f$
		StandardVariables operator+(const StandardVariables &b) const;

		/// @brief	Hashes monomial
		/// @return Hash of exponent vector (calls \ref generic_hasher)
		size_t operator()() const;
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief				Variables \f$e_1,...,e_n\f$ of degrees \f$|e_i|=i\f$
	///	@details			A monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ is stored as the vector \f$[a_1,...,a_n]\f$
	///	@tparam 	T 		The (integral) value type of the exponent vector.
	///	@tparam 	_deg 	The (integral) value type used in the degree function.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class T = int64_t, class _deg = int64_t>
	struct ElementarySymmetricVariables : public StandardVariables<T, _deg>
	{
		using StandardVariables<T, _deg>::StandardVariables;
		typedef _deg deg_t; ///<Degree typedef

		///	@brief	Computes degree of monomial on the \f$e_i\f$
		///	@return Degree \f$\sum_iia_i\f$ when ```*this```=\f$e_1^{a_1}\cdots e_n^{a_n}\f$
		deg_t degree() const;

		///	@brief		Returns the name of the variables \f$e_i\f$
		///	@return 	"e_i"
		///	@param 	i 	The variable index
		///	@param 	n 	The number variables
		static std::string name(int i, int n);
	};

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief					Factory class that provides the general interface of a generating basis for a subring of a polynomial ring.
	///	@details				Inherit from this class and construct data member \ref _generators through the child class (and optionally \ref generator_names and \ref generator_dimensions). \n
	///							The child class must also have a method \c find_exponent with signature:
	///							\code typename new_poly_t::exp_t find_exponent(const typename orig_poly_t::exp_t&); \endcode
	///							Example implementations are \c SymmetricBasis and \c TwistedChernBasis.
	///	@tparam spec_t  		Used for compile-time polymorphism (CRTP): must be the child class.
	///	@tparam orig_poly_t 	Type of polynomial on the original variables
	///	@tparam new_poly_t	 	Type of polynomial on the new variables (the \ref _generators)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class spec_t, class orig_poly_t, class new_poly_t>
	class PolynomialBasis
	{
	public:
		///	@brief		Transform a polynomial on the original variables to one on the generating basis
		/// @param 	a 	Polynomial on the original variables
		/// @return 	Polynomial on the new variables
		new_poly_t operator()(orig_poly_t a) const;

		///	@brief		Transform a polynomial on the generating basis into a polynomial on the original variables
		/// @param 	a 	Polynomial on the new variables
		/// @return 	Polynomial on the original variables
		orig_poly_t operator()(const new_poly_t &a) const;

		///	@brief		Constructor given number of variables
		/// @param num 	The number of variables for the polynomials
		PolynomialBasis(int num);

		///	@brief	Returns vector containing the generating basis
		const std::vector<orig_poly_t> &generators() const;

		///	@brief	Returns vector containing the dimensions of the generating basis (can be empty!)
		const std::vector<typename new_poly_t::deg_t> &dimensions() const;

		///	@brief	Returns vector containing the names of the generating basis (can be empty!)
		const std::vector<std::string> &names() const;

		///	@brief	The number of (the original) variables of the polynomial ring
		const int number_of_variables;

	protected:
		///	@brief	The generators of the polynomial basis, constructed in the inheriting class
		std::vector<orig_poly_t> _generators;

		///	@brief	The dimensions of the generators, optionally constructed in the inheriting class
		std::vector<typename new_poly_t::deg_t> generator_dimensions;

		///	@brief	The names of the generators, optionally constructed in the inheriting class
		std::vector<std::string> generator_names;

	private:
		orig_poly_t compute_product(const typename new_poly_t::exp_t &exponent) const;
	};

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief				Class for symmetric polynomials with no relations, allowing transformation from \f$x_i\f$ variables to \f$e_i\f$ variables and vice-versa.
	///	@tparam x_poly_t 	Type of Polynomial on the Standard_Variables \f$x_i\f$
	///	@tparam e_poly_t 	The of Polynomial on the ElementarySymmetricVariables \f$e_i\f$
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <class x_poly_t, class e_poly_t>
	class SymmetricBasis : public PolynomialBasis<SymmetricBasis<x_poly_t, e_poly_t>, x_poly_t, e_poly_t>
	{
	public:
		///	@brief		Constructor given number of variables
		///	@param num 	The number of variables for our symmetric polynomials
		SymmetricBasis(int num);

	private:
		using PolynomialBasis<SymmetricBasis<x_poly_t, e_poly_t>, x_poly_t, e_poly_t>::number_of_variables;
		using PolynomialBasis<SymmetricBasis<x_poly_t, e_poly_t>, x_poly_t, e_poly_t>::_generators;
		typedef typename x_poly_t::exp_t x_t;
		typedef typename e_poly_t::exp_t e_t;
		x_poly_t get_elementary_symmetric(int i) const;
		auto find_exponent(const x_t &term) const -> e_t;

		///	@brief	Befriending parent for CRTP.
		friend class PolynomialBasis<SymmetricBasis<x_poly_t, e_poly_t>, x_poly_t, e_poly_t>;
	};
}
#include "impl/Symmetric_Basis.ipp"