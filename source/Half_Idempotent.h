#pragma once
#include "Polynomials.h"
#include "Generators.h"

///@file
///@brief Contains the methods and classes for solving the following problem: If $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ produce minimal algebra generators for the fixed points of \f$R\f$ under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately), give an algorithm for writing a fixed point in terms of the generators
/// and an algorithm for producing the relations of those generators. 

namespace Symmetric_Polynomials {

	///Class for variables \f$x_1,...,x_n,y_1,...,y_n\f$ where \f$y_i^2=y_i\f$
	//
	///Parameter N is the number of variables in compile-time; set to -1 if unknown (default). Otherwise n=N.
	template<typename T, int N=-1>
	struct Half_Idempotent_Variables : public std::array<T,N> {
		using std::array<T,N>::array;
		///Constructor using size (does nothing as size is known in compile time for N!=-1)
		Half_Idempotent_Variables<T, N>(int) : std::array<T, N>() {}
		///The degree of a monomial on the \f$x_i,y_i\f$ with \f$|x_i|=1\f$ and \f$|y_i|=0\f$
		int degree() const {
			int degree = 0;
			for (int i = 0; i < N / 2; i++)
				degree += this->operator[](i);
			return degree;
		}
		///Multiplies monomials by adding their exponents, taking into account that \f$y_i^2=y_i\f$
		Half_Idempotent_Variables<T,N> operator+ (const Half_Idempotent_Variables<T,N>& other) const {
			Half_Idempotent_Variables<T,N> v;
			for (int i = 0; i < N; i++)
				v[i] = (*this)[i] + other[i];
			for (int i = N / 2; i < N; i++)
				v[i] = (v[i] > 0);
			return v;
		}
		///Divides monomials a/b with b|a by subtracting their exponents, taking into account that \f$y_i^2=y_i\f$
		Half_Idempotent_Variables<T,N> operator- (const Half_Idempotent_Variables<T, N>& other) const {
			Half_Idempotent_Variables<T, N> v;
			for (int i = 0; i < N; i++)
				v[i] = (*this)[i] - other[i];
			for (int i = N / 2; i < N; i++)
				v[i] = (v[i] > 0);
			return v;
		}
	};

	///Class specifying variables \f$x_1,...,x_n,y_1,...,y_n\f$ where \f$y_i^2=y_i\f$.
	//
	///The number of variables n is unknown in compile time.
	template<typename T>
	struct Half_Idempotent_Variables<T,-1>: public std::vector<T> {
		using std::vector<T>::vector;
		///The degree of a monomial on the x_i,y_i with |x_i|=1 and |y_i|=0
		int degree() const {
			int degree = 0;
			for (decltype(this->size()) i = 0; i < this->size() / 2; i++)
				degree += this->operator[](i);
			return degree;
		}
		///Multiplies monomials by adding their exponents, taking into account that y_i^2=y_i
		Half_Idempotent_Variables<T,-1> operator+ (const Half_Idempotent_Variables<T, -1>& other) const {
			Half_Idempotent_Variables<T, -1> v(this->size());
			for (size_t i = 0; i < v.size(); i++)
				v[i] = (*this)[i] + other[i];
			for (auto i = v.size() / 2; i < v.size(); i++)
				v[i] = (v[i] > 0);
			return v;
		}
		///Divides monomials a/b with b|a by subtracting their exponents, taking into account that y_i^2=y_i
		Half_Idempotent_Variables<T, -1> operator- (const Half_Idempotent_Variables<T, -1>& other) const {
			Half_Idempotent_Variables<T, -1> v(this->size());
			for (size_t i = 0; i < v.size(); i++)
				v[i] = (*this)[i] - other[i];
			for (auto i = v.size() / 2; i < v.size(); i++)
				v[i] = (v[i] > 0);
			return v;
		}
	};

	///Class for half-idempotent symmetric polynomials, allowing transformation from \f$x_i,y_i\f$ variables to \f$\alpha,c_i,\gamma_{s,i}\f$ variables and vice-versa.
	template<typename scalar_t, typename exponent_t>
	class Half_Idempotent_Basis {

		typedef std::vector<typename exponent_t::value_type> twistedChern_t;

	public:
		///Constructs the generators and the relation set
		Half_Idempotent_Basis(int n) : n(n) {
			set_generators();
			set_relations();
		}

		///Transform a symmetric polynomial on the \f$x_i,y_i\f$ into a polynomial on \f$\alpha,c_i,\gamma_{s,j}\f$
		Polynomial<scalar_t, twistedChern_t,1,1> operator()(Polynomial<scalar_t, exponent_t> a) const {
			Polynomial<scalar_t, twistedChern_t, 1, 1> decomposition(variable_names,dimensions);
			while (true) {
				twistedChern_t pwr_here(number_of_generators);
				auto max = a.highest_term();
				find_exponent(max.exponent(), pwr_here);
				auto product = this->compute_product(pwr_here);
				auto coeff = max.coeff() / product.highest_term().coeff();
				decomposition.insert(pwr_here, coeff);
				product *= coeff;
				if (a != product)
					a -= product;
				else
					return decomposition;
			}
		}

		///Transform a polynomial on the \f$\alpha,c_i,\gamma_{s,j}\f$ into a symmetric polynomial on the original variables \f$x_i,y_i\f$
		template<bool t, bool w>
		Polynomial<scalar_t, exponent_t> operator()(Polynomial<scalar_t, twistedChern_t,t,w> a) const {
			Polynomial<scalar_t, exponent_t> p(n);
			for (auto it = a.begin(); it != a.end(); ++it) {
				auto prod = compute_product(it.exponent());
				prod *= it.coeff();
				p += prod;
			}
			return p;
		}

		///The left hand side of the relations \f$\alpha,c_i,\gamma_{s,j}\f$ satisfy i.e. \f$\alpha^{n+1}, \alpha^s\gamma_{s,i}\f$ and \f$\gamma_{s,i}\gamma_{t,j}\f$ for \f$s<=t<=s+i\f$
		std::vector<Polynomial<scalar_t, twistedChern_t, 1, 1>> relations;

		Polynomial<scalar_t, exponent_t> Idem; ///<The generator \f$\alpha=y_1+...+y_n\f$
		std::vector<Polynomial<scalar_t, exponent_t>> Chern; ///<Vector consisting of the Chern classes \f$c_i\f$
		std::vector<std::vector<Polynomial<scalar_t, exponent_t>>> TwistedChern; ///<Vector consisting of the twisted Chern classes \f$\gamma_{s,j}\f$

	private:

		const int n; ///<Number of variables
		int number_of_generators; ///<Number of generators
		std::map<std::pair<int, int>, int> TwistedChern_indexer; ///<Transforms the usual pair index of a TwistedChern class to the single integer index we use for the vector TwistedChern


		///Computes a product of a monomial on our generators Idem, Chern, TwistedChern
		Polynomial<scalar_t, exponent_t> compute_product(const twistedChern_t& exponent) const {
			Polynomial<scalar_t, exponent_t> product = Idem ^ exponent[0];
			int counter = 1;
			for (const auto& c : Chern) {
				if (exponent[counter] != 0)
					product *= c ^ exponent[counter];
				counter++;
			}
			for (const auto& s : TwistedChern)
				for (const auto& i : s) {
					if (exponent[counter] != 0)
						product *= i ^ exponent[counter];
					counter++;
				}
			return product;
		}
		std::vector<int> dimensions;
		std::vector<std::string> variable_names; ///<The names of the generators Idem, Chern and TwistedChern

		Polynomial<scalar_t, exponent_t> get_Idem() {
			Polynomial<scalar_t, exponent_t> idem(2 * n);
			exponent_t mono(2 * n);
			for (int i = 2 * n - 1; i >= n; i--) {
				mono[i]=1;
				idem.insert(mono, 1);
				mono[i]=0;
			}
			return idem;
		}

		Polynomial<scalar_t, exponent_t> get_Chern(int i) {
			Polynomial<scalar_t, exponent_t> chern(2 * n);
			exponent_t mono(2 * n);
			auto c = Combination_Generator<typename exponent_t::value_type>(n, i);
			for (const auto& comb : c) {
				for (const auto& j : comb)
					mono[j] = 1;
				chern.insert(mono, 1);
				for (const auto& j : comb)
					mono[j] = 0;
			}
			return chern;
		}


		Polynomial<scalar_t, exponent_t> get_TwistedChern(int s, int i) {
			Polynomial<scalar_t, exponent_t> tchern(2 * n);
			exponent_t mono(2 * n);
			auto c_x = Combination_Generator<typename exponent_t::value_type>(n, s);
			auto c_y = Combination_Generator<typename exponent_t::value_type>(n - s, i);
			for (const auto& comb_x : c_x) {
				for (const auto& j : comb_x)
					mono[j]=1;
				std::vector<typename exponent_t::value_type> letters_not_in_comb_x;
				letters_not_in_comb_x.reserve(n - s);
				for (int j = 0; j < n; j++) {
					if (std::find(comb_x.begin(), comb_x.end(), j) == comb_x.end())
						letters_not_in_comb_x.push_back(j);
				}
				for (const auto& comb_y : c_y) {
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]]=1;
					tchern.insert(mono, 1);
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]]=0;
				}
				for (const auto& j : comb_x)
					mono[j] = 0;
			}
			return tchern;
		}


		void set_generators() {
			number_of_generators = 1 + (n * n + n) / 2;
			variable_names.reserve(number_of_generators);
			dimensions.reserve(number_of_generators);
			Idem = get_Idem();
			variable_names.push_back("a");
			dimensions.push_back(0);
			Chern.reserve(n);
			for (int i = 1; i <= n;i++) {
				Chern.push_back(get_Chern(i));
				variable_names.push_back("c_" + std::to_string(i));
				dimensions.push_back(i);
			}
			int index = 1 + n; //Idem,e_1,...,e_n
			TwistedChern.reserve((n * n - n) / 2);
			for (int s = 1; s <= n; s++) {
				std::vector<Polynomial<scalar_t, exponent_t>> cs;
				cs.reserve(n - s);
				for (int i = 1; i <= n - s; i++) {
					cs.push_back(get_TwistedChern(s, i));
					variable_names.push_back("c_{" + std::to_string(s) + "," + std::to_string(i) + "}");
					dimensions.push_back(s);
					TwistedChern_indexer[std::make_pair(s, i)] = index;
					index++;
				}
				TwistedChern.push_back(cs);
			}
		}

		void set_relations() {
			twistedChern_t comb(number_of_generators);
			comb[0] = n + 1;
			relations.emplace_back(comb,1, variable_names, dimensions); //Idem^{n+1} is a relation
			for (int s = 1; s < n; s++) {
				for (int i = 1; i <= n - s; i++) {
					std::fill(comb.begin(), comb.end(), 0);
					comb[0] = s;
					comb[TwistedChern_indexer.at(std::make_pair(s, i))] = 1; // Idem^s TwistedChern_{s,i} is a relation
					relations.emplace_back(comb,1,variable_names,dimensions);
				}
			}
			for (int s = 1; s < n; s++)
				for (int t = s; t < n; t++)
					for (int i = std::max(t - s, 1); i <= n - s; i++)
						for (int j = 1; j <= n - t; j++) {
							if (s == t && j < i)
								continue;
							std::fill(comb.begin(), comb.end(), 0);
							comb[TwistedChern_indexer.at(std::make_pair(s, i))]++;
							comb[TwistedChern_indexer.at(std::make_pair(t, j))]++; // TwistedChern_{s,i} TwistedChern_{t,j} is is a relation if s<=t<=s+i
							relations.emplace_back(comb, 1, variable_names, dimensions);
						}
		}

		void find_exponent(const exponent_t& term, twistedChern_t& exponent) const {
			if (term[n] > 0) {//clear the consecutive Idem's at the start
				int consecutive_u_at_start = 0;
				exponent_t us_at_start(2 * n);
				for (int i = n; i < 2*n; i++) {
					if (term[i] > 0) {
						us_at_start[i]=1;
						consecutive_u_at_start++;
					}
					else
						break;
				}
				exponent[0] = consecutive_u_at_start;
				find_exponent(term - us_at_start, exponent);
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
			if (last_consecutive_u_from_end == 2 * n) {//no Idem's
				for (int i = 1; i <= n; i++) {
					if (i == n)
						exponent[i] += term[i - 1];
					else
						exponent[i] += term[i - 1] - term[i];
				}
				return;
			}
			if (last_consecutive_u_from_end > n) {//not first Idem
				//add c_{s,i} to the exponent
				exponent[TwistedChern_indexer.at(std::pair(last_consecutive_u_from_end - n,number_consecutive_u))]++;
				find_exponent(term - TwistedChern[last_consecutive_u_from_end - n - 1][number_consecutive_u - 1].highest_term().exponent(), exponent);
			}
		}

	};

	///Prints all relations in the description of the fixed points of R=Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i) using Idem, Chern and TwistedChern 
	template<typename scalar_t, typename exponent_t>
	void print_half_idempotent_relations(int n, bool print =0, bool verify = 0) {
		Half_Idempotent_Basis<scalar_t, exponent_t> hib(n);
//#pragma omp parallel for num_threads(12)
		for (const auto & rel:hib.relations) {
			auto p= hib(rel);
			auto q = hib(p);
			if (print)
				std::cout << rel << " = " << q << "\n";
			if (verify)
				std::cout << "Verification: " << (p== hib(q)) << "\n";
		}
	}
}
