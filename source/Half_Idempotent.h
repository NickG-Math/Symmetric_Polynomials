#pragma once
#include "Permutations.h"
#include <map>

///@file
///@brief If $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$, the goal is to produce algebra generators for the fixed points of R under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately), give an algorithm for writing a fixed point in terms of the generators
/// and an algorithm for producing the relations of those generators

namespace Symmetric_Polynomials {
	////////////////////////////////////////////////////////////////////////////////////////////
	///The class containing the generators for the fixed points of $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ under the \f$\Sigma_n\f$ permutation action.
	//
	/// There are three types of generators (the following is the default order): Idem (\f$y_1+...+y_n\f$), Chern_1,...,Chern_n (elementary symmetric on x_1,...,x_n) as well as the TwistedChern_{s,j} indexed on 1<=s<=n and 1<=j<=n-s
	/// The latter involve both the x_i and y_i. The order on the TwistedChern is lexicographic on s,j.
	///////////////////////////////////////////////////////////////////////////////////
	template<typename scalar_t>
	class half_idempotent_generators {

	public:
		///Constructs the generators and the relation set
		half_idempotent_generators(int n) : n(n) {
			set_generators();
			set_relations();
		}
		const int n; ///<Number of variables
		int number_of_generators; ///<Number of generators
		std::map<std::pair<int, int>, int> TwistedChern_indexer; ///<Transforms the usual pair index of a TwistedChern class to the single integer index we use for the vector TwistedChern

		///Computes a product of a monomial on our generators Idem, Chern, TwistedChern
		polynomial<scalar_t, halfidempotent> compute_product(const std::vector<int>& exponent) const {
			polynomial<scalar_t, halfidempotent> product = power(Idem, exponent[0]);
			int counter = 1;
			for (int i = 0; i < Chern.size(); i++) {
				if (exponent[counter] != 0)
					product = multiply<scalar_t, halfidempotent>(product, power(Chern[i], exponent[counter]));
				counter++;
			}
			for (const auto& s : TwistedChern)
				for (const auto& i : s) {
					if (exponent[counter] != 0)
						product = multiply(product, power(i, exponent[counter]));
					counter++;
				}
			return product;
		}

		///Print a polynomial on the generators Idem, Chern and TwistedChern
		std::string print(const polynomial<scalar_t, relations_base>& m) const {
			return m.print(variable_names);
		}

		std::vector<std::string> variable_names; ///<The names of the generators Idem, Chern and TwistedChern
		std::vector<std::vector<int>> relations; ///<The left hand side of the relations Idem, Chern and TwistedChern satisfy i.e. Idem^{n+1}, Idem^sTwistedChern_{s,i} and TwistedChern_{s,i}TwistedChern_{t,j} for s<=t<=s+i

		polynomial<scalar_t, halfidempotent> Idem; ///<a=y_1+...+y_n
		std::vector<polynomial<scalar_t, halfidempotent>> Chern; ///<Vector consisting of the Chern classes c_i
		std::vector<std::vector<polynomial<scalar_t, halfidempotent>>> TwistedChern; ///<Vector consisting of the twisted Chern classes b_{s,j}.

	private:
		polynomial<scalar_t, halfidempotent> get_Idem() {
			std::vector< monomial<scalar_t, halfidempotent>> Idem_v;
			Idem_v.reserve(n);
			std::vector<int> mono(2 * n);
			for (int i = 2 * n - 1; i >= n; i--) {
				mono[i] = 1;
				Idem_v.push_back(monomial<scalar_t, halfidempotent>(1, mono));
				mono[i] = 0;
			}
			return polynomial<scalar_t, halfidempotent>(Idem_v);
		}

		polynomial<scalar_t, halfidempotent> get_Chern(int i) {
			std::vector< monomial<scalar_t, halfidempotent>> Chern_v;
			Chern_v.reserve(binomial(n, i));
			std::vector<int> mono(2 * n);
			auto c = combinations_generator(n, i);
			for (const auto& comb : c) {
				for (const auto& j : comb)
					mono[j] = 1;
				Chern_v.push_back(monomial<scalar_t, halfidempotent>(1, mono));
				for (const auto& j : comb)
					mono[j] = 0;
			}
			return polynomial<scalar_t, halfidempotent>(Chern_v.rbegin(), Chern_v.rend(), 1);
		}


		polynomial<scalar_t, halfidempotent> get_TwistedChern(int s, int i) {
			std::vector< monomial<scalar_t, halfidempotent>> TChern_v;
			TChern_v.reserve(binomial(n, s) * binomial(n - s, i));
			std::vector<int> mono(2 * n);
			auto c_x = combinations_generator(n, s);
			auto c_y = combinations_generator(n - s, i);
			for (const auto& comb_x : c_x) {
				for (const auto& j : comb_x)
					mono[j] = 1;
				std::vector<char> letters_not_in_comb_x;
				letters_not_in_comb_x.reserve(n - s);
				for (int j = 0; j < n; j++) {
					if (std::find(comb_x.begin(), comb_x.end(), j) == comb_x.end())
						letters_not_in_comb_x.push_back(j);
				}
				for (const auto& comb_y : c_y) {
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]] = 1;
					TChern_v.push_back(monomial<scalar_t, halfidempotent>(1, mono));
					for (const auto& j : comb_y)
						mono[n + letters_not_in_comb_x[j]] = 0;
				}
				for (const auto& j : comb_x)
					mono[j] = 0;
			}
			return polynomial<scalar_t, halfidempotent>(TChern_v.rbegin(), TChern_v.rend(), 1);
		}


		void set_generators() {
			number_of_generators = 1 + (n * n + n) / 2;
			variable_names.reserve(number_of_generators);
			Idem=get_Idem();
			variable_names.push_back("a");
			Chern.reserve(n);
			for (int i = 1; i <= n;i++) {
				Chern.push_back(get_Chern(i));
				variable_names.push_back("c_" + std::to_string(i));
			}
			int index = 1 + n; //Idem,e_1,...,e_n
			TwistedChern.reserve((n * n - n) / 2);
			for (int s = 1; s <= n; s++) {
				std::vector<polynomial<scalar_t, halfidempotent>> cs;
				cs.reserve(n - s);
				for (int i = 1; i <= n - s; i++) {
					cs.push_back(get_TwistedChern(s, i));
					variable_names.push_back("b_{" + std::to_string(s) + "," + std::to_string(i) + "}");
					TwistedChern_indexer[std::make_pair(s, i)] = index;
					index++;
				}
				TwistedChern.push_back(cs);
			}
		}

		void set_relations() {
			std::vector<int> comb(number_of_generators);
			comb[0] = n + 1;
			relations.push_back(comb); //Idem^{n+1} is a relation
			for (int s = 1; s <= n; s++) {
				for (int i = 1; i <= n - s; i++) {
					std::fill(comb.begin(), comb.end(), 0);
					comb[0] = s;
					comb[TwistedChern_indexer.at(std::make_pair(s, i))] = 1; // Idem^s TwistedChern_{s,i} is a relation
					relations.push_back(comb);
				}
			}
			for (int s = 1; s <= n; s++)
				for (int t = s; t <= n; t++)
					for (int i = std::max(t - s,1); i <= n - s; i++)
						for (int j = 1; j <= n - t; j++) {
							if (s==t && j<i)
								continue;
							std::fill(comb.begin(), comb.end(), 0);
							comb[TwistedChern_indexer.at(std::make_pair(s, i))]++;
							comb[TwistedChern_indexer.at(std::make_pair(t, j))]++; // TwistedChern_{s,i} TwistedChern_{t,j} is is a relation if s<=t<=s+i
							relations.push_back(comb);
						}
		}
	};


	///Class for writing every element in the fixed points of $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ as a polynomial on the generators Idem, Chern and TwistedChern (a "decomposition")
	template<typename scalar_t>
	struct decomposition_half_idempotent : public half_idempotent_generators<scalar_t> {
		polynomial<scalar_t, relations_base> decomposition; ///< The polynomial on the Idem, Chern and TwistedChern generators. We use relations_base since the relations for these generators are very complicated

		///Constructor given number of variables n
		decomposition_half_idempotent(int n) : half_idempotent_generators<scalar_t>(n) {}

		///Initialize and decompose into polynomial on Idem,Chern and TwistedChern
		decomposition_half_idempotent(const polynomial<scalar_t, halfidempotent>& a) : decomposition_half_idempotent<scalar_t>(a.monos[0].exponent.size() / 2) { decompose(a); }


		///Computes the decomposition (polynomial on the generators Idem, Chern and TwistedChern)
		void decompose(const polynomial<scalar_t, halfidempotent>& a) {
			decomposition.monos.clear();
			decompose_recursive(a);
		}

		///Verifies that the decomposition is correct (i.e. expanding the generators in the original variables \f$x_i,y_i\f$ gives the original polynomial).
		bool verify(const polynomial<scalar_t, halfidempotent>& a) const {
			auto result = decomposition.monos[0].coeff * this->compute_product(decomposition.monos[0].exponent);
			for (int i = 1; i < decomposition.monos.size(); i++)
				result = result + decomposition.monos[i].coeff * this->compute_product(decomposition.monos[i].exponent);
			return (a == result);
		}

	private:


		void find_exponent(const monomial<scalar_t, halfidempotent>& term, std::vector<int>& exponent) {
			if (term.exponent[this->n] > 0) {//clear the consecutive Idem's at the start
				int consecutive_u_at_start = 0;
				std::vector<int> us_at_start(2 * this->n);
				for (int i = this->n; i < term.exponent.size(); i++) {
					if (term.exponent[i] > 0) {
						us_at_start[i] = 1;
						consecutive_u_at_start++;
					}
					else
						break;
				}
				exponent[0] = consecutive_u_at_start;
				find_exponent(term / monomial<scalar_t, halfidempotent>(1, us_at_start), exponent);
				return;
			}
			bool found_last_u = 0;
			int last_consecutive_u_from_end = 2 * this->n;
			int number_consecutive_u = 0;
			for (int i = 2 * this->n - 1; i >= this->n; i--) {
				if (term.exponent[i] == 0 && found_last_u)
					break;
				else if (term.exponent[i] > 0) {
					last_consecutive_u_from_end = i;
					number_consecutive_u++;
					found_last_u = 1;
				}
			}
			if (last_consecutive_u_from_end == 2 * this->n) {//no Idem's
				for (int i = 1; i <= this->n; i++) {
					if (i == this->n)
						exponent[i] += term.exponent[i - 1];
					else
						exponent[i] += term.exponent[i - 1] - term.exponent[i];
				}
				return;
			}
			if (last_consecutive_u_from_end > this->n) {//not first Idem
				//add c_{s,i} to the exponent
				exponent[this->TwistedChern_indexer[std::make_pair(last_consecutive_u_from_end - this->n, number_consecutive_u)]]++;
				find_exponent(term / this->TwistedChern[last_consecutive_u_from_end - this->n - 1][number_consecutive_u - 1].highest_term(), exponent);
			}
		}

		void decompose_recursive(const polynomial<scalar_t, halfidempotent>& a) {
			std::vector<int> pwr_here(this->number_of_generators);
			auto max = a.highest_term();
			find_exponent(max, pwr_here);
			auto product = this->compute_product(pwr_here);
			auto coeff = max.coeff / product.highest_term().coeff;
			decomposition.monos.push_back(monomial<scalar_t, relations_base>(coeff, pwr_here));
			auto coeff_product = coeff * product;
			if (a != coeff_product)
				decompose_recursive(a - coeff_product);
		}

	};

	///Prints polynomial
	template<typename scalar_t>
	std::ostream& operator<<(std::ostream& os, const decomposition_half_idempotent<scalar_t>& dec) {
		os << dec.print(dec.decomposition);
		return os;
	}

	///Prints all relations in the description of the fixed points of R=Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i) using Idem, Chern and TwistedChern 
	template<typename scalar_t>
	void print_half_idempotent_relations(int n, bool verify=0) {
		decomposition_half_idempotent<scalar_t> dec(n);
		for (const auto& rel : dec.relations) {
			auto polynomial_relation = dec.compute_product(rel);
			dec.decompose(polynomial_relation);
			std::cout << dec.print(monomial<scalar_t, relations_base>(1, rel)) << " = " << dec << "\n"; 
			if (verify)
				std::cout << "Verification: " << dec.verify(polynomial_relation) << "\n";
		}
	}
}
