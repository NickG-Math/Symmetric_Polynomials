#pragma once
#include "Symmetric_Basis.hpp"

///@file
///@brief Contains the methods and classes for solving the following problem: If $$R=\mathbb Z[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ produce minimal algebra generators for the fixed points of \f$R\f$ under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately), give an algorithm for writing a fixed point in terms of the generators
/// and an algorithm for producing the relations of those generators. 


namespace Symmetric_Polynomials {

	///Wrapping array and vector in the same interface
	//
	///@tparam T The value type of the array/vector
	///@tparam N The size if it's an array or 0 if it's a vector
	template<typename T, size_t N = 0>
	struct Array_Vector_Wrapper;

	///The variables \f$x_1,...,x_n,y_1,...,y_n\f$ where \f$y_i^2=y_i\f$ and \f$|x_i|=1\f$, \f$|y_i|=0\f$
	//
	///Monomial \f$x_1^{a_1}\cdots x_n^{a_n}y_1^{a_{n+1}}\cdots y_n^{a_{2n}}\f$ is stored as vector/array \f$[a_1,...,a_{2n}]\f$
	///@tparam T 	The (integral) value type of the exponent vector.
	///@tparam _deg The (integral) value type used in the degree function.
	///@tparam N 	The number of variables in compile-time; set to 0 if unknown (default). Otherwise N=\f$2n\f$.
	template<typename T = int64_t, typename _deg = int64_t, size_t N = 0>
	struct Half_Idempotent_Variables : public Array_Vector_Wrapper<T, N> {
		using Array_Vector_Wrapper<T, N>::Array_Vector_Wrapper;

		///Multiplies monomials by adding their exponents
		//
		///@return \f$[a_1+b_1,...,a_n+b_n, \max(a_{n+1},b_{n+1}), ..., \max(a_{2n},b_{2n})]\f$ where *this=\f$[a_1,...,a_{2n}]\f$
		///@param b Second exponent \f$[b_1,...,b_{2n}]\f$
		Half_Idempotent_Variables operator+ (const Half_Idempotent_Variables&) const;

		///Divides monomials by subtracting their exponents
		//
		///@return \f$[a_1-b_1,...,a_n-b_n, |a_{n+1}-b_{n+1}|, ..., |a_{2n}-b_{2n}|]\f$ where *this=\f$[a_1,...,a_{2n}]\f$
		///@param b Second exponent \f$[b_1,...,b_{2n}]\f$. We must have \f$b_i\le a_i\f$ for every \f$i\f$.
		Half_Idempotent_Variables operator- (const Half_Idempotent_Variables&) const;

		typedef _deg deg_t;					///<Degree typedef

		///Computes degree of monomial on on the \f$x_i,y_i\f$ with \f$|x_i|=1\f$ and \f$|y_i|=0\f$ 
		//
		///@return \f$\sum_{i=1}^na_i\f$ for monomial \f$x_1^{a_1}\cdots y_n^{a_{2n}}\f$ (*this=\f$[a_1,...,a_{2n}]\f$)
		deg_t degree() const;				

		///Returns the names of the variables \f$x_i,y_i\f$
		//
		///@return "x_i" if i<n and "y_{i-n}" if i>n
		///@param i The variable index
		///@param num The number variables = \f$2n\f$									
		static std::string name(int i, int num);	

		///Returns hash of monomial
		size_t operator ()() const;
	};

	
	///The twisted Chern generators as variables \f$\gamma_{s,j}\f$
	//
	///Monomial \f$\prod_{s,j}\gamma_{s,j}^{a_{s_j}}\f$ is stored as vector \f$[a_{0,1},...,a_{0,n},a_{1,0},a_{1,1},...,a_{n-1,1},a_{n,0}]\f$
	///
	///This class does NOT provide functions for degrees or variable names: these are provided as pointers directly in Half_Idempotent_Basis
	///@tparam T 	The (integral) value type of the exponent vector.
	///@tparam _deg The (integral) value type used in the degree function.
	template<typename T = int64_t, typename _deg = int64_t>
	struct Twisted_Chern_Variables : public std::vector<T> {
		using std::vector<T>::vector;
		typedef _deg deg_t;														  ///<Degree typedef

		///<Multiplies monomials by adding their exponents. 
		///
		///No relations are used as these wouldn't produce monomials
		///@return \f$[a_{0,1}+b_{0,1},...,a_{n,0}+b_{n,0}]\f$ where *this=\f$[a_{0,1},...,a_{n,0}]\f$
		///@param b Second exponent \f$[b_{0,1},...,b_{n,0}]\f$
		Twisted_Chern_Variables operator+ (const Twisted_Chern_Variables& b) const; 
		size_t operator ()() const;												  ///<Returns hash of monomial
	};

	///Class for half-idempotent symmetric polynomials, allowing transformation from \f$x_i,y_i\f$ variables to \f$\gamma_{s,i}\f$ variables and vice-versa.
	//
	///@tparam _xy_container The container type on the Half_Idempotent_Variables \f$x_i,y_i\f$
	///@tparam _chern_container The container type on the Twisted_Chern_Variables \f$\gamma_{s,j}\f$
	template<typename _xy_container, typename _chern_container>
	class Half_Idempotent_Basis : public Polynomial_Basis<Half_Idempotent_Basis<_xy_container,_chern_container>, _xy_container, _chern_container> {
	public:
		///Constructs the generators and the relation set given \f$n\f$ in \f$x_1,...,x_n,y_1,...,y_n\f$.
		//
		///@param n \f$n\f$ is half(!) the number of variables
		Half_Idempotent_Basis(int n);

		///Returns const& of vector containing the Polynomials \f$\alpha_s\gamma_{s,i}\f$ and \f$\gamma_{s,i}\gamma_{t,j}\f$ for \f$0<s<=t<=s+i\f$ and \f$i,j>0\f$ (the LHS of the relations satisfied by the generators)
		const auto& relations() const;

		///Returns generator \f$\gamma_{s,j}\f$.
		//
		/// @return \f$\gamma_{s,j}\f$ as Polynomial on the \f$x_i,y_i\f$ variables
 		/// @param s The index of \f$s\f$ of \f$\gamma_{s,j}\f$.
 		/// @param j The index of \f$j\f$ of \f$\gamma_{s,j}\f$.
		const auto& generator(int s, int j) const;

		typedef _xy_container xy_container_t; ///<The container of polynomials on x,y variables
		typedef _chern_container chern_container_t; ///<The container of polynomials on x,y variables

		///Befriending parent for CRTP.
		friend class Polynomial_Basis<Half_Idempotent_Basis<_xy_container, _chern_container>, _xy_container, _chern_container>;

	private:

		typedef typename Polynomial<_xy_container>::exp_t xy_t;
		typedef typename Polynomial<_chern_container>::exp_t chern_t;
		using Polynomial_Basis<Half_Idempotent_Basis<_xy_container, _chern_container>, _xy_container, _chern_container>::_generators;
		using Polynomial_Basis<Half_Idempotent_Basis<_xy_container, _chern_container>, _xy_container, _chern_container>::generator_names;
		using Polynomial_Basis<Half_Idempotent_Basis<_xy_container, _chern_container>, _xy_container, _chern_container>::generator_dimensions;

		const int n; //Half the variable number
		const int number_of_generators;  //Number of c_{s,j}
		std::vector<Polynomial<_chern_container>> _relations;

		std::map<std::array<int,2>, int> generator_double_index;  //Takes (s,j) to the index in the _generators vector
		int index(int s, int j) const; //Transforms index (s,j) of c_{s,j} into the corresponding index in the generators vector.
		auto create_generator(int s, int i);
		void set_generators();
		void set_relations();
		auto find_exponent(const xy_t& term) const;
		void find_exponent_recursive(const xy_t& term, chern_t& exponent) const;
	};

	///Aliases Half_Idempotent_Basis for certain default template parameters
	///@tparam scl_t The type of scalars eg int
	///@tparam exp_value_t The value type of the exponent vectors
	///@tparam deg_t The type of the degree of the monomials
	///@tparam ordered Whether to use std::map or std::unordered_map
	template<typename scl_t, typename exp_value_t = int64_t, typename deg_t = int64_t, bool ordered = 1>
	using Half_Idempotent_Basis_Default = Half_Idempotent_Basis<default_container<scl_t, Half_Idempotent_Variables<exp_value_t, deg_t>, ordered>, default_container<scl_t, Twisted_Chern_Variables<exp_value_t, deg_t>, ordered>>;

	///Prints all relations in the description of the fixed points of \f$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)\f$ in terms of \f$\alpha_i, c_i, \gamma_{s,j}\f$ (printed as a_i,c_i,c_{s,j} in the console) 
	//
	///@tparam HIB The Half_Idempotent_Basis type
	///@param n Half the number of variables (the \f$n\f$)
	///@param print Whether we want to print the relations to the console
	///@param verify Whether to verify the relations
	///@param verify_verbose Whether to verify and print the verification to the console
	template<typename HIB = Half_Idempotent_Basis_Default<int64_t, uint64_t, uint64_t, 0>>
	void print_half_idempotent_relations(int n, bool print = 0, bool verify = 0, bool verify_verbose = 0);

}
#include "impl/Half_Idempotent.ipp"