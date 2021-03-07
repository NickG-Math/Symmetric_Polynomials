#pragma once
#include "Symmetric_Basis.h"

///@file
///@brief Contains the methods and classes for solving the following problem: If $$R=\mathbb Z[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ produce minimal algebra generators for the fixed points of \f$R\f$ under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately), give an algorithm for writing a fixed point in terms of the generators
/// and an algorithm for producing the relations of those generators. 


namespace {
	//wrapping array and vector in the same interface

	template<typename T, size_t N = 0>
	struct array_vector_wrapper : public std::array<T, N>{
		using std::array<T, N>::array;
		array_vector_wrapper(size_t) : std::array<T, N>() {}
	};

	template<typename T>
	struct array_vector_wrapper<T, 0> : public std::vector<T> {
		using std::vector<T>::vector;
	};

}

namespace Symmetric_Polynomials {


	///Class for variables \f$x_1,...,x_n,y_1,...,y_n\f$ where \f$y_i^2=y_i\f$
	//
	///Parameter N is the number of variables in compile-time; set to -1 (i.e. uint64_t max) if unknown (default). Otherwise n=N.
	template<typename T = int64_t, typename _deg = int64_t, size_t N = 0>
	struct Half_Idempotent_Variables : public array_vector_wrapper<T, N> {
		using array_vector_wrapper<T, N>::array_vector_wrapper;
		typedef _deg deg_t;///<Degree typedef

		///The degree of a monomial on the \f$x_i,y_i\f$ with \f$|x_i|=1\f$ and \f$|y_i|=0\f$
		deg_t degree() const {
			deg_t degree = 0;
			for (size_t i = 0; i < this->size() / 2; i++)
				degree += this->operator[](i);
			return degree;
		}
		///Multiplies monomials by adding their exponents, taking into account that \f$y_i^2=y_i\f$
		Half_Idempotent_Variables operator+ (const Half_Idempotent_Variables& other) const {
			Half_Idempotent_Variables v(this->size());
			for (size_t i = 0; i < this->size(); i++)
				v[i] = (*this)[i] + other[i];
			for (size_t i = this->size() / 2; i < this->size(); i++)
				v[i] = (v[i] > 0);
			return v;
		}
		///Divides monomials a/b when b|a by subtracting their exponents, taking into account that \f$y_i^2=y_i\f$
		Half_Idempotent_Variables operator- (const Half_Idempotent_Variables& other) const {
			Half_Idempotent_Variables v(this->size());
			for (size_t i = 0; i < this->size(); i++)
				v[i] = (*this)[i] - other[i];
			for (size_t i = this->size() / 2; i < this->size(); i++)
				v[i] = (v[i] > 0);
			return v;
		}
		///Names of variables x_i,y_i
		static std::string name(int i, int number_of_variables) {
			if (i < number_of_variables / 2)
				return "x_" + std::to_string(i + 1);
			else
				return "y_" + std::to_string(i - number_of_variables / 2 + 1);
		}
		///Standard hash function
		size_t operator ()() const {
			return generic_hasher(*this);
		}
	};

	///Class specifying the generators \f$\alpha_i,c_i,\gamma_{s,j}\f$
	//
	///Does not include degrees or variable names:  these are provided as pointers in Half_Idempotent_Basis
	template<typename T = int64_t, typename _deg = int64_t>
	struct Twisted_Chern_Variables : public std::vector<T> {
		using std::vector<T>::vector;
		typedef _deg deg_t;///<Degree typedef
		Twisted_Chern_Variables operator+ (const Twisted_Chern_Variables& other) const {
			Standard_Variables v(*this);
			for (size_t i = 0; i < this->size(); i++)
				v[i] += other[i];
			return v;
		}
		///Standard hash function
		size_t operator ()() const {
			return generic_hasher(*this);
		}
	};

	///Class for half-idempotent symmetric polynomials, allowing transformation from \f$x_i,y_i\f$ variables to \f$\gamma_{s,i}\f$ variables and vice-versa.
	template<typename, typename>
	class Half_Idempotent_Basis;

	///Class for half-idempotent symmetric polynomials, allowing transformation from \f$x_i,y_i\f$ variables to \f$\gamma_{s,i}\f$ variables and vice-versa.
	template<typename xy_container_t, typename chern_container_t>
	class Half_Idempotent_Basis : public Polynomial_Basis<Half_Idempotent_Basis<xy_container_t,chern_container_t>, xy_container_t, chern_container_t> {

		typedef typename Polynomial<xy_container_t>::exp_t xy_t;
		typedef typename Polynomial<chern_container_t>::exp_t chern_t;

		using Polynomial_Basis<Half_Idempotent_Basis<xy_container_t,chern_container_t>, xy_container_t, chern_container_t>::_generators;
		using Polynomial_Basis<Half_Idempotent_Basis<xy_container_t,chern_container_t>, xy_container_t, chern_container_t>::generator_names;
		using Polynomial_Basis<Half_Idempotent_Basis<xy_container_t,chern_container_t>, xy_container_t, chern_container_t>::generator_dimensions;

	public:

		///Constructs the generators and the relation set given \f$n\f$ in \f$x_1,...,x_n,y_1,...,y_n\f$. Warning: \f$n\f$ is half(!) the number of variables
		Half_Idempotent_Basis(int n) : Polynomial_Basis<Half_Idempotent_Basis<xy_container_t,chern_container_t>, xy_container_t, chern_container_t> (2*n), n(n), number_of_generators(n + (n * n + n) / 2) {
			set_generators();
			set_relations();
		}

		///Returns const& of vector containing the Polynomials \f$\alpha_s\gamma_{s,i}\f$ and \f$\gamma_{s,i}\gamma_{t,j}\f$ for \f$0<s<=t<=s+i\f$ and \f$i,j>0\f$ (the LHS of the relations satisfied by the generators)
		const auto& relations() const {
			return _relations;
		}

		///Returns const& of \f$\gamma_{s,j}\f$ as Polynomial on the \f$x_i,y_i\f$ variables. 
		const auto& generator(int s, int j) const {
			return _generators[index(s,j)];
		}

	private:

		const int n; //Half the variable number
		const int number_of_generators;  //Number of c_{s,j}

		std::vector<Polynomial<chern_container_t>> _relations;

		std::map<std::array<int,2>, int> generator_double_index;  //Takes (s,j) to the index in the _generators vector
		//Transforms index \f$(s,j)\f$ of \f$\gamma_{s,j}\f$ into the corresponding index in the generators vector.
		int index(int s, int j) const {
			return generator_double_index.at({ s,j });
		}


		auto create_generator(int s, int i) {
			Polynomial<xy_container_t> tchern(2 * n);
			xy_t mono(2 * n);
			auto c_x = Combination_Generator<typename xy_t::value_type>(n, s);
			auto c_y = Combination_Generator<typename xy_t::value_type>(n - s, i);
			for (const auto& comb_x : c_x) {
				for (const auto& j : comb_x)
					mono[j] = 1;
				std::vector<typename xy_t::value_type> letters_not_in_comb_x;
				letters_not_in_comb_x.reserve(n - s);
				for (int j = 0; j < n; j++) {
					if (std::find(comb_x.begin(), comb_x.end(), j) == comb_x.end())
						letters_not_in_comb_x.push_back(j);
				}
				for (const auto& comb_y : c_y) {
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]] = 1;
					tchern.insert(mono, 1);
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]] = 0;
				}
				for (const auto& j : comb_x)
					mono[j] = 0;
			}
			return tchern;
		}


		void set_generators() {
			generator_names.reserve(number_of_generators);
			generator_dimensions.reserve(number_of_generators);
			_generators.reserve(number_of_generators);
			int ind = 0;
			for (int s = 0; s <= n; s++) {
				for (int i = 0; i <= n - s; i++) {
					if (i == 0 && s == 0)
						continue;
					_generators.push_back(create_generator(s, i));
					if (i == 0)
						generator_names.push_back("c_" + std::to_string(s));
					else if (s== 0)
						generator_names.push_back("a_" + std::to_string(i));
					else
						generator_names.push_back("c_{" + std::to_string(s) + "," + std::to_string(i) + "}");
					generator_dimensions.push_back(s);
					generator_double_index[{s, i}] = ind;
					ind++;
				}
			}
		}

		void set_relations() {
			for (int s = 0; s < n; s++)
				for (int t = s; t < n; t++)
					for (int i = 1; i <= n - s; i++)
						for (int j = 1; j <= n - t; j++) {
							if (t > s + i || (s == t && j < i) || (s==0 && t!=i) ) //if t>s+i no relation, if s==t && j<i we have a symmetric relation, if s==0 && t!=i we can get if from the relation between a_j and a_1
								continue;
							chern_t comb(number_of_generators);
							comb[generator_double_index[{s, i}]]++;
							comb[generator_double_index[{t, j}]]++; 
							_relations.emplace_back(comb, 1, generator_dimensions.data(), generator_names.data());
						}
		}

		chern_t find_exponent(const xy_t& term) const {
			chern_t exponent(number_of_generators);
			find_exponent_recursive(term, exponent);
			return exponent;
		}

		void find_exponent_recursive(const xy_t& term, chern_t& exponent) const {
			if (term[n] > 0) {//clear the consecutive y_i at the start
				int consecutive_u_at_start = 0;
				xy_t us_at_start(2 * n);
				for (int i = n; i < 2 * n && term[i]>0; i++) {
					us_at_start[i] = 1;
					consecutive_u_at_start++;
				}
				if (consecutive_u_at_start > 0)
					exponent[index(0,consecutive_u_at_start)] = 1;
				find_exponent_recursive(term - us_at_start, exponent);
				return;
			}
			bool found_last_u = 0;
			int last_consecutive_u_from_end = 2 * n;
			int number_consecutive_u = 0;
			for (int i = 2 * n - 1; i >= n; i--) {
				if (term[i] == 0 && found_last_u)
					break;
				else if (term[i] > 0) {
					last_consecutive_u_from_end = i;
					number_consecutive_u++;
					found_last_u = 1;
				}
			}
			if (last_consecutive_u_from_end == 2 * n) {//no y_i means it's an elementary symmetric on the x_i
				for (int i = 1; i < n; i++) 
					exponent[index(i, 0)] += term[i - 1] - term[i];
				exponent[index(n, 0)] += term[n - 1];
				return;
			}
			if (last_consecutive_u_from_end > n) {//not first y_i so add c_{s,i} to the exponent
				auto ind = index(last_consecutive_u_from_end - n, number_consecutive_u);
				exponent[ind]++;
				find_exponent_recursive(term - _generators[ind].highest_term().exponent(), exponent);
			}
		}

		///Befriending parent for CRTP.
		friend class Polynomial_Basis<Half_Idempotent_Basis<xy_container_t, chern_container_t>, xy_container_t, chern_container_t>;

	};

	///Aliases Half_Idempotent_Basis for certain default template parameters
	template<typename scl_t, typename exp_value_t = int64_t, typename deg_t = int64_t, bool ordered = 1>
	using Half_Idempotent_Basis_Default = Half_Idempotent_Basis<default_container<scl_t, Half_Idempotent_Variables<exp_value_t, deg_t>, ordered>, default_container<scl_t, Twisted_Chern_Variables<exp_value_t, deg_t>, ordered>>;

	///Prints all relations in the description of the fixed points of \f$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)\f$ in terms of \f$\alpha_i, c_i, \gamma_{s,j}\f$ (printed as a_i,c_i,c_{s,j} in the console) 
	template<typename HIB = Half_Idempotent_Basis_Default<int64_t, uint64_t, uint64_t, 0>>
	void print_half_idempotent_relations(int n, bool print = 0, bool verify = 0, bool verify_verbose=0) {
		HIB hib(n);
		//#pragma omp parallel for num_threads(12)
		for (const auto& rel : hib.relations()) {
			auto p = hib(rel);
			auto q = hib(p);
			if (print) 
				std::cout << rel << " = " << q << "\n";
			if (verify) {
				if (p != hib(q))
					throw("Verification failed!");
				else {
					std::cout << "Verified!";
					if (verify_verbose)
						std::cout << "In x, y variables both LHS and RHS are : " << p;
					std::cout << "\n";
				}
			}
		}
	}

}
