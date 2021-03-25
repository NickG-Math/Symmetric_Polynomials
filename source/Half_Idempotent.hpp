#pragma once
#include "Symmetric_Basis.hpp"
#include <sstream>

/////////////////////////////////////////////////////////////////////////
///	@file
///	@brief 		Contains the methods and classes for symmetric polynomials with half idempotent variables
///	@details	The goal is to solve the following problem:
///				If \f[R=\mathbf Z[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)\f]
///				produce minimal algebra generators for the fixed points of \f$R\f$
///				under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately),
/// 			give an algorithm for writing a fixed point in terms of the generators
///				and an algorithm for producing the relations of those generators.
/////////////////////////////////////////////////////////////////////////

namespace symmp
{

	/////////////////////////////////////////////////////////////////////////
	///	@brief		Wrapping array and vector in the same interface
	///	@tparam T	The value type of the array/vector
	///	@tparam N	The size if it's an array or 0 if it's a vector
	/////////////////////////////////////////////////////////////////////////
	template <typename T, size_t N = 0>
	struct ArrayVectorWrapper;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief			Variables \f$x_1,...,x_n,y_1,...,y_n\f$ with \f$y_i^2=y_i\f$ and \f$|x_i|=1\f$, \f$|y_i|=0\f$
	///	@details		Monomial \f$x_1^{a_1}\cdots x_n^{a_n}y_1^{a_{n+1}}\cdots y_n^{a_{2n}}\f$ is stored as vector/array \f$[a_1,...,a_{2n}]\f$
	///	@tparam T 		The (integral) value type of the exponent vector.
	///	@tparam _deg	The (integral) value type used in the degree function.
	///	@tparam N 		The number of variables in compile-time; set to 0 if unknown (default). Otherwise N=\f$2n\f$.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename T = int64_t, typename _deg = int64_t, size_t N = 0>
	struct HalfIdempotentVariables : public ArrayVectorWrapper<T, N>
	{
		using ArrayVectorWrapper<T, N>::ArrayVectorWrapper;

		/// @brief		Multiplies monomials by adding their exponents
		/// @param  b	Second exponent \f$[b_1,...,b_{2n}]\f$
		/// @return		Exponent \f$[a_1+b_1,...,a_n+b_n, \max(a_{n+1},b_{n+1}), ..., \max(a_{2n},b_{2n})]\f$ where *this=\f$[a_1,...,a_{2n}]\f$
		HalfIdempotentVariables operator+(const HalfIdempotentVariables &b) const;

		///	@brief		Divides monomials by subtracting their exponents
		///	@param	b	Second exponent \f$[b_1,...,b_{2n}]\f$. We must have \f$b_i\le a_i\f$ for every \f$i\f$.
		///	@return		Exponent \f$[a_1-b_1,...,a_n-b_n, |a_{n+1}-b_{n+1}|, ..., |a_{2n}-b_{2n}|]\f$ where *this=\f$[a_1,...,a_{2n}]\f$
		HalfIdempotentVariables operator-(const HalfIdempotentVariables &b) const;

		typedef _deg deg_t; ///<Degree typedef

		///	@brief	Computes degree of monomial on on the \f$x_i,y_i\f$ with \f$|x_i|=1\f$ and \f$|y_i|=0\f$
		///	@return	Degree \f$\sum_{i=1}^na_i\f$ for monomial \f$x_1^{a_1}\cdots y_n^{a_{2n}}\f$ (*this=\f$[a_1,...,a_{2n}]\f$)
		deg_t degree() const;

		///	@brief			Returns the names of the variables \f$x_i,y_i\f$
		///	@param i		The variable index
		///	@param num		The number variables = \f$2n\f$
		///	@return			"x_i" if i<n and "y_{i-n}" if i>n
		static std::string name(int i, int num);

		/// @brief	Hashes monomial
		///	@return	Hash of exponent vector (calls generic_hasher)
		size_t operator()() const;
	};

	///	@brief			The twisted Chern generators as variables \f$\gamma_{s,j}\f$
	///	@details		Monomial \f$\prod_{s,j}\gamma_{s,j}^{a_{s_j}}\f$ is stored as vector \f$[a_{0,1},...,a_{0,n},a_{1,0},a_{1,1},...,a_{n-1,1},a_{n,0}]\f$
	///	@note			This class does NOT provide functions for degrees or variable names: these are provided as pointers directly in TwistedChernBasis
	///	@tparam T 		The (integral) value type of the exponent vector.
	///	@tparam _deg	The (integral) value type used in the degree function.
	template <typename T = int64_t, typename _deg = int64_t>
	struct TwistedChernVariables : public std::vector<T>
	{
		using std::vector<T>::vector;
		typedef _deg deg_t; ///<Degree typedef

		///	@brief		Multiplies monomials by adding their exponents.
		///	@note		No relations are used as these wouldn't produce monomials
		///	@param b	Second exponent \f$[b_{0,1},...,b_{n,0}]\f$
		///	@return		Exponent \f$[a_{0,1}+b_{0,1},...,a_{n,0}+b_{n,0}]\f$ where *this=\f$[a_{0,1},...,a_{n,0}]\f$
		TwistedChernVariables operator+(const TwistedChernVariables &b) const;

		/// @brief	Hashes monomial
		/// @return Hash of exponent vector (calls generic_hasher)
		size_t operator()() const;
	};

	///	@brief						Class for half-idempotent symmetric polynomials, allowing transformation from \f$x_i,y_i\f$ variables to \f$\gamma_{s,i}\f$ variables and vice-versa.
	///	@tparam _xy_poly_t			The container type on the HalfIdempotentVariables \f$x_i,y_i\f$
	///	@tparam _chern_poly_t		The container type on the TwistedChernVariables \f$\gamma_{s,j}\f$
	template <typename _xy_poly_t, typename _chern_poly_t>
	class TwistedChernBasis : public PolynomialBasis<TwistedChernBasis<_xy_poly_t, _chern_poly_t>, _xy_poly_t, _chern_poly_t>
	{
	public:
		///	@brief		Constructs the generators and the relation set given \f$n\f$ in \f$x_1,...,x_n,y_1,...,y_n\f$.
		///	@param n	Half(!) the number of variables
		///	@attention	The parameter \c n is half(!) the number of variables (the \f$n\f$ in \f$BU(n)\f$)
		TwistedChernBasis(int n);

		///	@brief		Stores the relations \f$\gamma_{s,i}\gamma_{t,j}\f$ for \f$0<s<=t<=s+i\f$ and \f$i,j>0\f$
		///	@return		const& of vector containing the relation Polynomials \f$\gamma_{s,i}\gamma_{t,j}\f$
		const auto &relations() const;

		///	@brief		Returns generator \f$\gamma_{s,j}\f$.
		/// @param s	The index of \f$s\f$ of \f$\gamma_{s,j}\f$.
		/// @param j	The index of \f$j\f$ of \f$\gamma_{s,j}\f$.
		/// @return		The polynomial \f$\gamma_{s,j}\f$ on the \f$x_i,y_i\f$ variables
		const auto &generator(int s, int j) const;

		///	@brief	Befriending parent for CRTP.
		friend class PolynomialBasis<TwistedChernBasis<_xy_poly_t, _chern_poly_t>, _xy_poly_t, _chern_poly_t>; 

	private:
		typedef typename _xy_poly_t::exp_t xy_t;
		typedef typename _chern_poly_t::exp_t chern_t;
		using PolynomialBasis<TwistedChernBasis<_xy_poly_t, _chern_poly_t>, _xy_poly_t, _chern_poly_t>::_generators;
		using PolynomialBasis<TwistedChernBasis<_xy_poly_t, _chern_poly_t>, _xy_poly_t, _chern_poly_t>::generator_names;
		using PolynomialBasis<TwistedChernBasis<_xy_poly_t, _chern_poly_t>, _xy_poly_t, _chern_poly_t>::generator_dimensions;

		const int n;					//Half the variable number
		const int number_of_generators; //Number of c_{s,j}
		std::vector<_chern_poly_t> _relations;

		std::map<std::array<int, 2>, int> generator_double_index; //Takes (s,j) to the index in the _generators vector
		int index(int s, int j) const;							  //Transforms index (s,j) of c_{s,j} into the corresponding index in the generators vector.
		auto create_generator(int s, int i);
		void set_generators();
		void set_relations();
		auto find_exponent(const xy_t &term) const;
		void find_exponent_recursive(const xy_t &term, chern_t &exponent) const;
	};

	///	@brief					Prints all relations in the description of the fixed points of \f$R=\mathbf Z[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)\f$ in terms of \f$\alpha_i, c_i, \gamma_{s,j}\f$ (printed as ``` a_i, c_i, c_{s,j} ``` in the console)
	///	@tparam xy_poly_t		The type of polynomial on the \f$x_i,y_i\f$ variables
	///	@tparam chern_poly_t	The type of polynomial on the \f$\gamma_{s,j}\f$ variables
	///	@param n				Half the number of variables (the \f$n\f$)
	///	@param print			Whether we want to print the relations to the console
	///	@param verify			Whether to verify the relations
	///	@param verify_verbose	Whether to verify and print the verification to the console
	template <typename xy_poly_t, typename chern_poly_t>
	void print_half_idempotent_relations(int n, bool print = 0, bool verify = 0, bool verify_verbose = 0);

}
#include "impl/Half_Idempotent.ipp"
