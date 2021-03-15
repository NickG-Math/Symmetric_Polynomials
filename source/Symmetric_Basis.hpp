#pragma once
#include "Polynomials.hpp"
#include "Generators.hpp"

///@file
///@brief Contains methods for writing symmetric polynomials (with no relations) in terms of the elementary symmetric ones. Also contains a factory class that generalizes the interface to any subring of a polynomial ring

namespace Symmetric_Polynomials {

	///The standard variables \f$x_i\f$ in a polynomial, with \f$|x_i|=1\f$ and no relations.
	//
	///A monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ is stored as the vector \f$[a_1,...,a_n]\f$
	///@tparam T 	The (integral) value type of the exponent vector.
	///@tparam _deg The (integral) value type used in the degree function.
	template<typename T = int64_t, typename _deg = int64_t>
	struct Standard_Variables : public std::vector<T> {
		using std::vector<T>::vector;
		typedef _deg deg_t;												///<Degree typedef

		///Computes degree of monomial on standard variables \f$x_i\f$
		//
		///@return \f$\sum_ia_i\f$ for monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ (*this=\f$[a_1,...,a_n]\f$)
		deg_t degree() const;	

		///Returns the names of the standard variables \f$x_i\f$
		//
		///@return "x_i"
		///@param i The variable index
		///@param n The number variables									
		static std::string name(int i, int n);							

		///Multiplies monomials by adding their exponents
		//
		///@return \f$[a_1+b_1,...,a_n+b_n]\f$ where *this=\f$[a_1,...,a_n]\f$
		///@param b Second exponent \f$[b_1,...,b_n]\f$
		Standard_Variables operator+ (const Standard_Variables& b) const; 

		///Returns hash of monomial
		size_t operator ()() const;
	};

	///Variables \f$e_1,...,e_n\f$ denoting the elementary symmetric polynomials \f$e_i=\sigma_i\f$ of degrees \f$|e_i|=i\f$
	//
	///A monomial \f$x_1^{a_1}\cdots x_n^{a_n}\f$ is stored as the vector \f$[a_1,...,a_n]\f$
	///@tparam T 	The (integral) value type of the exponent vector.
	///@tparam _deg The (integral) value type used in the degree function.
	template<typename T = int64_t, typename _deg = int64_t>
	struct Elementary_Symmetric_Variables : public Standard_Variables<T, _deg> {
		using Standard_Variables<T, _deg>::Standard_Variables;
		typedef _deg deg_t;					///<Degree typedef
		
		///Computes degree of monomial on the \f$e_i\f$
		//
		///@return \f$\sum_iia_i\f$ for monomial \f$e_1^{a_1}\cdots e_n^{a_n}\f$
		deg_t degree() const;	

		///Returns the name of the variables \f$e_i\f$
		//
		///@return "e_i"
		///@param i The variable index
		///@param n The number variables												
		static std::string name(int, int);
	};

	///////////////////////////////////////////////////////////////////////////////////
	///Factory class that provides the general interface of a generating basis for a subring of a polynomial ring. Example implementations are Symmetric_Basis and Half_Idempotent_Basis.
	//
	///Inherit from this class and construct data member _generators through the child class (and optionally generator_names and generator_dimensions).
	///
	///The child class must also have a method find_exponent with singature: Polynomial<container_2_t>::exp_t find_exponent(const Polynomial<container_1_t>::exp_t&);
	///@tparam specialization  Used for compile-time polymorphism (CRTP): set it to be the child class. 
	///@tparam container_1_t The container type of the original variables (of the polynomial ring)
	///@tparam container_2_t The container type of the new variables that are the basis elements (of the subring)
	///////////////////////////////////////////////////////////////////////////////////
	template<typename specialization, typename container_1_t, typename container_2_t>
	class Polynomial_Basis {
	public:
		///Transform a polynomial on the original variables to one on the generating basis
		//
		/// @return Polynomial on the new variables
 		/// @param a Polynomial on the original variables
		Polynomial<container_2_t> operator()(Polynomial<container_1_t> a) const;

		///Transform a polynomial on the generating basis into a polynomial on the original variables
		//
		/// @return Polynomial on the original variables
 		/// @param a Polynomial on the new variables
		Polynomial<container_1_t> operator()(const Polynomial<container_2_t>& a) const;

		///Constructor given number of variables
		//
		/// @param number_of_variables The number of variables for the polynomials
		Polynomial_Basis(int number_of_variables);

		///Returns vector containing the generating basis
		const auto& generators() const;

		///Returns vector containing the dimensions of the generating basis (can be empty!)
		const auto& dimensions() const;

		///Returns vector containing the names of the generating basis (can be empty!)
		const auto& names() const;

		///The number of (the original) variables of the polynomial ring
		const int number_of_variables;

	protected:
		///The generators of the polynomial basis, constructed in the inheriting class
		std::vector<Polynomial<container_1_t>> _generators;

		///The dimensions of the generators, optionally constructed in the inheriting class
		std::vector<typename Polynomial<container_2_t>::deg_t> generator_dimensions;

		///The names of the generators, optionally constructed in the inheriting class
		std::vector<std::string> generator_names;

	private:
		Polynomial<container_1_t> compute_product(const typename Polynomial<container_2_t>::exp_t& exponent) const;
	};

	///Class for symmetric polynomials with no relations, allowing transformation from \f$x_i\f$ variables to \f$e_i\f$ variables and vice-versa.
	//
	///@tparam x_container_t The container type on the Standard_Variables \f$x_i\f$
	///@tparam e_container_t The container type on the Elementary_Symmetric_Variables \f$e_i\f$
	template<typename x_container_t, typename e_container_t>
	class Symmetric_Basis : public Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t> {
	public:
		///Constructor given number of variables
		//
		///@param number_of_variables The number of variables for our symmetric polynomials
		Symmetric_Basis(int number_of_variables);
	private:
		using Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>, x_container_t, e_container_t>::number_of_variables;
		using Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>, x_container_t, e_container_t>::_generators;
		typedef typename Polynomial<x_container_t>::exp_t x_t;
		typedef typename Polynomial<e_container_t>::exp_t e_t;
		Polynomial<x_container_t> get_elementary_symmetric(int i) const;
		auto find_exponent(const x_t& term) const;

		///Befriending parent for CRTP.
		friend class Polynomial_Basis<Symmetric_Basis<x_container_t, e_container_t>,x_container_t,e_container_t>;
	};


	///Aliases Symmetric_Basis for certain default template parameters
	///@tparam scl_t The type of scalars eg int
	///@tparam exp_value_t The value type of the exponent vectors
	///@tparam deg_t The type of the degree of the monomials
	///@tparam ordered Whether to use std::map or std::unordered_map
	template<typename scl_t, typename exp_value_t = int64_t, typename deg_t = int64_t, bool ordered = 1>
	using Symmetric_Basis_Default = Symmetric_Basis<default_container<scl_t, Standard_Variables<exp_value_t, deg_t>, ordered>, default_container<scl_t, Elementary_Symmetric_Variables<exp_value_t, deg_t>, ordered>>;

}
#include "impl/Symmetric_Basis.ipp"