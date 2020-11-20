#pragma once
#include "Permutations.h"
#include <map>

///@file
///@brief If $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$, the goal is to produce algebra generators for the fixed points of R under the \f$\Sigma_n\f$ action (permuting the \f$x_i,y_i\f$ separately), give an algorithm for writing a fixed point in terms of the generators
/// and an algorithm for producing the relations of those generators

namespace Symmetric_Polynomials{
////////////////////////////////////////////////////////////////////////////////////////////
///The class containing the generators for the fixed points of $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ under the \f$\Sigma_n\f$ permutation action.
//
/// There are three types of generators (the following is the default order): Idem (\f$y_1+...+y_n\f$), Chern_1,...,Chern_n (elementary symmetric on x_1,...,x_n) as well as the TwistedChern_{s,j} indexed on 1<=s<=n and 1<=j<=n-s
/// The latter involve both the x_i and y_i. The order on the TwistedChern is lexicographic on s,j.
///////////////////////////////////////////////////////////////////////////////////
template<typename scalar_t>
class half_idempotent_generators {

public:
	///Constructs the generators, given the number of variables n
	half_idempotent_generators(int n) : half_idempotent_generators(n, all_e_indices(n), all_e_indices(n - 1), all_c_si_indices(n)) {}


	const int n; ///<Number of variables
	int number_of_generators; ///<Number of generators


	std::map<std::pair<int, int>, int> TwistedChern_indexer; ///<Transforms the usual pair index of a TwistedChern class to the single integer index we use for the vector TwistedChern


	std::vector<int> dimensions; ///< The dimensions of the generators Idem, Chern, TwistedChern in this order

	///Computes a product of a monomial on our generators Idem, Chern, TwistedChern
	polynomial<scalar_t, halfidempotent> compute_product(const std::vector<int>& pwr) const {
		polynomial<scalar_t, halfidempotent> product = power(polynomial<scalar_t, halfidempotent>(orbit(Idem, permutations)), pwr[0]);
		int counter = 1;
		for (int i = 0; i < Chern.size(); i++) {
			if (pwr[counter] != 0)
				product = multiply<scalar_t, halfidempotent>(product, power(polynomial<scalar_t, halfidempotent>(orbit(Chern[i], permutations)), pwr[counter]));
			counter++;
		}
		for (const auto& s : TwistedChern)
			for (const auto& i : s) {
				if (pwr[counter] != 0)
					product = multiply(product, power(polynomial<scalar_t, halfidempotent>(orbit(i, permutations)), pwr[counter]));
				counter++;
			}
		product.condense();
		return product;
	}

	///Computes the dominant term of a product of a monomial on our generators Idem, Chern, TwistedChern
	monomial<scalar_t, halfidempotent> compute_product_max_term(const std::vector<int>& pwr) const {
		//The product of Idem,TwistedChern is computed completely. 
		//The max term of the product of Idem,Chern,TwistedChern must then be the max term of Chern*max term of Idem,TwistedChern

		std::vector<int> product_e_max_vector(2 * n);
		int counter = 1;
		for (int i = 0; i < Chern.size(); i++) {
			if (pwr[counter] != 0)
				product_e_max_vector = product_e_max_vector + pwr[counter] * Chern[i].powers;
			counter++;
		}
		monomial<scalar_t, halfidempotent> product_e_max(1, product_e_max_vector);

		polynomial<scalar_t, halfidempotent> product_u_c = power(polynomial<scalar_t, halfidempotent>(orbit(Idem, permutations)), pwr[0]);
		for (const auto& s : TwistedChern)
			for (const auto& i : s) {
				if (pwr[counter] != 0)
					product_u_c = multiply(product_u_c, power(polynomial<scalar_t, halfidempotent>(orbit(i, permutations)), pwr[counter]));
				counter++;
			}
		auto product_max = product_e_max * product_u_c.highest_term();
		return product_max;
	}


	///Print a polynomial on the generators Idem, Chern and TwistedChern
	std::string print(const polynomial<scalar_t,relations_base>& m) const {
		return m.print(variable_names);
	}

protected:

	monomial<scalar_t, halfidempotent> Idem; ///The highest term of y_1+...+y_n
	std::vector<monomial<scalar_t, halfidempotent>> Chern; ///Vector consisting of the highest term of the Chern classes
	std::vector<std::vector<monomial<scalar_t, halfidempotent>>> TwistedChern; ///Vector consisting of the highest term of the twisted Chern classes.


	///Returns the max powers in a product of Idem, Chern and TwistedChern that has the desired degree
	std::vector<int> max_powers_given_degree(int degree) const {
		std::vector<int> max;
		max.reserve(number_of_generators);
		max.push_back(n);
		for (const auto i : e_indices)
			max.push_back(degree / i);
		for (int t = 0; t < c_s_indices.size(); t++) {
			auto s = c_s_indices[t];
			for (const auto i : c_si_indices[t]) {
				max.push_back(std::min(degree / s, 1));
			}
		}
		return max;
	}

	const std::vector<std::vector<char>> permutations; ///<All permutations in n variables


private:
	std::vector<std::string> variable_names; ///<The names of the generators Idem, Chern and TwistedChern
	std::vector<std::vector<char>> relations;

	const std::vector<int> e_indices, c_s_indices;
	const std::vector<std::vector<int>> c_si_indices;

	half_idempotent_generators(int n, const std::vector<int>& e_indices, const std::vector<int>& c_s_indices, const std::vector<std::vector<int>>& c_si_indices) :
		n(n), permutations(all_permutations(n)), e_indices(e_indices), c_s_indices(c_s_indices), c_si_indices(c_si_indices) {
		set_generators();
		set_relations();
	}

	static std::vector<int> all_e_indices(int n) {
		std::vector<int> e_indices(n);
		std::iota(e_indices.begin(), e_indices.end(), 1);
		return e_indices;
	}

	static std::vector<std::vector<int>> all_c_si_indices(int n) {
		std::vector<std::vector<int>> c_si_indices;
		c_si_indices.reserve(n);
		for (int s = 1; s < n; s++)
			c_si_indices.push_back(all_e_indices(n - s));
		return c_si_indices;
	}


	void set_generators() {
		number_of_generators = 1 + (n * n + n) / 2;
		dimensions.reserve(number_of_generators);
		variable_names.reserve(number_of_generators);
		Idem=monomial<scalar_t, halfidempotent>(1, std::vector<int>(2 * n));
		Idem.powers[n] = 1;
		dimensions.push_back(0);
		variable_names.push_back("a");
		Chern.reserve(e_indices.size());
		for (const auto i : e_indices) {
			std::vector<int> prod_of_first_i(2 * n);
			for (int j = 0; j < i; j++)
				prod_of_first_i[j] = 1;
			Chern.push_back(monomial<scalar_t, halfidempotent>(1, prod_of_first_i));
			dimensions.push_back(i);
			variable_names.push_back("c_" + std::to_string(i));
		}
		int index = 1 + e_indices.size(); //Idem,e_1,...,e_n
		TwistedChern.reserve(c_s_indices.size());
		for (int t = 0; t < c_s_indices.size(); t++) {
			auto s = c_s_indices[t];
			std::vector<monomial<scalar_t, halfidempotent>> cs;
			cs.reserve(c_si_indices[t].size());
			for (const auto i : c_si_indices[t]) {
				std::vector<int> csi(2 * n);
				for (int j = 0; j < s; j++)
					csi[j] = 1;
				for (int j = n + s; j < n + s + i; j++)
					csi[j] = 1;
				cs.push_back(monomial<scalar_t, halfidempotent>(1, csi));
				dimensions.push_back(s);
				variable_names.push_back("b_{" + std::to_string(s) + "," + std::to_string(i) + "}");
				TwistedChern_indexer[std::make_pair(s, i)] = index;
				index++;
			}
			TwistedChern.push_back(cs);
		}
	}

	void set_relations() {
		std::vector<char> comb(number_of_generators);
		for (int t = 0; t < c_s_indices.size(); t++) {
			auto s = c_s_indices[t];
			for (const auto i : c_si_indices[t]) {
				std::fill(comb.begin(), comb.end(), 0);
				comb[0] = s;
				comb[TwistedChern_indexer[std::make_pair(s, i)]] = 1; // Idem^s TwistedChern_{s,i} is a relation
				relations.push_back(comb);
			}
		}
		for (int ind = 0; ind < c_s_indices.size(); ind++)
			for (int r = ind; r < c_s_indices.size(); r++)
				for (const auto i : c_si_indices[ind]){
					auto s = c_s_indices[ind];
					auto t = c_s_indices[r];
					if (t <= s + i) {
						for (const auto j : c_si_indices[r]) {
							std::fill(comb.begin(), comb.end(), 0);
							if (s != t || i != j) { //TwistedChern^2 is not accounted for by the max vector
								comb[TwistedChern_indexer[std::make_pair(s, i)]] = 1;
								comb[TwistedChern_indexer[std::make_pair(t, j)]] = 1; // TwistedChern_{s,i} TwistedChern_{t,j} is is a relation if s<=t<=s+i
								relations.push_back(comb);
							}
						}
					}
				}
	}
};


///Class for writing every element in the fixed points of $$R=\mathbb Q[x_1,...,x_n,y_1,...,y_n]/(y_i^2=y_i)$$ as a polynomial on the generators Idem, Chern and TwistedChern (a "decomposition")
template<typename scalar_t>
struct decomposition_half_idempotent : public half_idempotent_generators<scalar_t> {
	polynomial<scalar_t,relations_base> decomposition; ///< The polynomial on the Idem, Chern and TwistedChern generators. We use relations_base since the relations for these generators are very complicated

	///Constructor given number of variables n
	decomposition_half_idempotent(int n) : half_idempotent_generators<scalar_t>(n) {}

	///Initialize and decompose into polynomial on Idem,Chern and TwistedChern
	decomposition_half_idempotent(const polynomial<scalar_t, halfidempotent>& a) : decomposition_half_idempotent<scalar_t>(a.monos[0].powers.size()/2) { decompose(a); }


	///Computes the decomposition (polynomial on the generators Idem, Chern and TwistedChern)
	void decompose(const polynomial<scalar_t, halfidempotent>& a) {
		decomposition.monos.clear();
		decompose_recursive(a);
	}

	///Verifies that the decomposition is correct (i.e. expanding the generators in the original variables \f$x_i,y_i\f$ gives the original polynomial).
	bool verify(const polynomial<scalar_t, halfidempotent>& a) const {
		auto result = decomposition.monos[0].coeff * this->compute_product(decomposition.monos[0].powers);
		for (int i = 1; i < decomposition.monos.size(); i++)
			result = result + decomposition.monos[i].coeff * this->compute_product(decomposition.monos[i].powers);
		return (a == result);
	}

private:


	void find_power(const monomial<scalar_t, halfidempotent>& term, std::vector<int>& power) {
		if (term.powers[this->n] > 0) {//clear the consecutive Idem's at the start
			int consecutive_u_at_start = 0;
			std::vector<int> us_at_start(2 * this->n);
			for (int i = this->n; i < term.powers.size(); i++) {
				if (term.powers[i] > 0) {
					us_at_start[i] = 1;
					consecutive_u_at_start++;
				}
				else
					break;
			}
			power[0] = consecutive_u_at_start;
			find_power(term / monomial<scalar_t, halfidempotent>(1, us_at_start), power);
			return;
		}
		bool found_last_u = 0;
		int last_consecutive_u_from_end = 2 * this->n;
		int number_consecutive_u = 0;
		for (int i = 2 * this->n - 1; i >= this->n; i--) {
			if (term.powers[i] == 0 && found_last_u)
				break;
			else if (term.powers[i] > 0) {
				last_consecutive_u_from_end = i;
				number_consecutive_u++;
				found_last_u = 1;
			}
		}
		if (last_consecutive_u_from_end == 2 * this->n) {//no Idem's
			for (int i = 1; i <= this->n; i++) {
				if (i == this->n)
					power[i] += term.powers[i - 1];
				else
					power[i] += term.powers[i - 1] - term.powers[i];
			}
			return;
		}
		if (last_consecutive_u_from_end > this->n) {//not first Idem
			//add c_{s,i} to the power
			power[this->TwistedChern_indexer[std::make_pair(last_consecutive_u_from_end - this->n, number_consecutive_u)]]++;
			find_power(term / this->TwistedChern[last_consecutive_u_from_end - this->n - 1][number_consecutive_u - 1], power);
		}
	}

	void decompose_recursive(const polynomial<scalar_t, halfidempotent>& a) {
		std::vector<int> pwr_here(this->dimensions.size());
		auto max = max_in_orbit(a.highest_term(), this->permutations);
		find_power(max, pwr_here);
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
void print_half_idempotent_relations(int n) {
	decomposition_half_idempotent<scalar_t> dec(n);
	std::vector<int> pwr(dec.dimensions.size());
	pwr[0] = n + 1;
	auto a = dec.compute_product(pwr);
	dec.decompose(a);
	std::cout << dec.print(monomial<scalar_t,relations_base>(1,pwr)) << " = " << dec << "\n"; //<< " Verified: " << dec.verify(a) << "\n";

	for (int s = 0; s < n-1; s++) {
		for (int i = 0; i < n-s-1; i++) {
			for (int t = s; t < std::min(s + i, n-1); t++)
				for (int j = 0; j < n-t-1; j++) {
					std::fill(pwr.begin(), pwr.end(), 0);
					pwr[dec.TwistedChern_indexer[std::make_pair(s + 1, i + 1)]]++;
					pwr[dec.TwistedChern_indexer[std::make_pair(t + 1, j + 1)]]++;
					auto a = dec.compute_product(pwr);
					dec.decompose(a);
					std::cout << dec.print(monomial<scalar_t, relations_base>(1, pwr)) << " = " << dec << "\n"; //" Verified: " << dec.verify(a) << "\n";
				}
		}
	}
}
}
