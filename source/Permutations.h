#pragma once
#include "Generators.h"

///@file
///@brief Contains methods for producing permutations, orbits of monomials under symmetric group action, and writing symmetric polynomials in terms of elementary symmetric ones

namespace Symmetric_Polynomials {

	///Returns $\binom{n,m}$
	long binomial(int top, int bot) {
		//safe implementation
		if (bot > top - bot)
			bot = top - bot;
		int binom = top;
		for (int i = 1; i <= bot-1; i++) { // $n/1 * (n-1)/2 * \cdots (n-k+1)/k$
			binom *= (top - i);
			binom /= i+1;
		}
		return binom;
	}


	///Returns $n!$
	long factorial(int n) {
		if (n == 0 || n == 1)
			return 1;
		else
			return n * factorial(n - 1);
	}

	////////////////////////////////////////////////////////////////////////////////////////
	///Generates all permutations on a number of letters
	//
	//Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a permutations_generator object. Then i will be a permutation
	//////////////////////////////////////////////////////////////////////////////////
	class permutations_generator {
		const char n;
	public:
		///Returns the total nuumber of permutations i.e. the factorial
		long size() const {
			return factorial(n);
		}
		///Constructor using number of letters
		permutations_generator(char n) : n(n) {}
		///Constant iterator that is used in a ranged for loop to generate the permutations. Non constant version is illegal
		class const_iterator : public generating_const_iterator<permutations_generator::const_iterator, std::vector<char>> {
			void update() {
				if (!std::next_permutation(generated.begin(), generated.end()))
					completed = 1;
			}
			const_iterator() {}
			friend class permutations_generator;
			friend class generating_const_iterator<permutations_generator::const_iterator, std::vector<char>>;
		};
		///Initial generator
		const_iterator begin() const {
			const_iterator it;
			it.generated.resize(n);
			std::iota(it.generated.begin(), it.generated.end(), 0);
			it.completed = 0;
			return it;
		}
		///Terminal generator.
		const_iterator end() const {
			return generating_const_iterator<permutations_generator::const_iterator, std::vector<char>>::terminal();
		}

	};



	////////////////////////////////////////////////////////////////////////////////////////
///Generates all combinations on a number of letters making a number of choices
//
//Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a combinations_generator object. Then i will be a combination
//////////////////////////////////////////////////////////////////////////////////
	class combinations_generator {
		const int total, choices;
	public:
		///Total number of combinations
		long size() const {
			return binomial(total,choices);
		}
		///Constructor using the total number of elements and the amount of choices we make
		combinations_generator(int total, int choices) : total(total), choices(choices) {}

		///Constant iterator that is used in a ranged for loop to generate the combinations. Non constant version is illegal
		class const_iterator : public generating_const_iterator<combinations_generator::const_iterator, std::vector<char>> {
			void update() {
				for (int i = choices - 1; i >= 0; i--) {
					if (generated[i] - i < total - choices) {
						generated[i]++;
						for (int j = i + 1; j < choices; j++) {
							generated[j] = generated[j - 1] + 1;
						}
						return;
					}
				}
				completed = 1;
				return;
			}
			const_iterator() {}
			friend class combinations_generator;
			friend class generating_const_iterator<combinations_generator::const_iterator, std::vector<char>>;
			int total, choices;
		};
		///Initial generator
		const_iterator begin() const{
			const_iterator it;
			it.generated.resize(choices);
			std::iota(it.generated.begin(), it.generated.end(), 0);
			it.total = total;
			it.choices = choices;
			it.completed=0;
			return it;
		}
		///Terminal generator.
		const_iterator end() const {
			return generating_const_iterator<combinations_generator::const_iterator, std::vector<char>>::terminal();
		}
	};


	///Returns vector of all permutations on n letters
	std::vector<std::vector<char>> all_permutations(char n) {
		std::vector<std::vector<char>> v;
		v.reserve(factorial(n));
		const permutations_generator p(n);
		for (const auto& i : p)
			v.push_back(i);
		return v;
	}

	///Returns vector of all combinations on n letters choosing m many  
	std::vector<std::vector<char>> all_permutations(char n, char m) {
		std::vector<std::vector<char>> v;
		v.reserve(binomial(n,m));
		const combinations_generator p(n,m);
		for (const auto& i : p)
			v.push_back(i);
		return v;
	}


	///Computes the orbit of a monomial under a set of permutations
	template<typename scalar_t, typename rel_t, typename perm_t>
	std::vector<monomial<scalar_t, rel_t>> orbit(const monomial<scalar_t, rel_t>& a, const perm_t& perms) {
		return orbit(a, perms, perms.size());
	}

	///Computes the orbit of a monomial under a set of permutations. The number_of_orbits is provided as an argument to limit memory usage
	template<typename scalar_t, typename rel_t, typename perm_t>
	std::vector<monomial<scalar_t, rel_t>> orbit(const monomial<scalar_t, rel_t>& a, const perm_t& perms, int number_of_orbits) {
		std::vector<monomial<scalar_t, rel_t>> o;
		o.reserve(number_of_orbits);
		for (const auto& i : perms) {
			auto permed = permute(a, i);
			auto it = std::find(o.begin(), o.end(), permed);
			if (it == o.end())
				o.push_back(permed);
		}
		return o;
	}

	///Returns the greatest monomial in an orbit of permutations
	template<typename scalar_t, typename rel_t, typename perm_t>
	monomial<scalar_t, rel_t> max_in_orbit(const monomial<scalar_t, rel_t>& a, const perm_t& perms) {
		auto max = a;
		for (const auto& i : perms) {
			auto permed = permute(a, i);
			if (max < permed)
				max = permed;
		}
		return max;
	}

	///Returns a vector of all monomials in given variables and degree organized in orbits under the symmetric group action
	template<typename scalar_t, typename rel_t>
	std::vector<std::vector<monomial<scalar_t, rel_t>>> all_monomial_orbits(int variables, int degree) {
		auto monos = monomial_basis<scalar_t, rel_t>(variables, degree);
		std::vector<std::vector<monomial<scalar_t, rel_t>>> orbits;
		std::vector<long> already_done;
		orbits.reserve(monos.size());
		already_done.reserve(monos.size());
		orbits.reserve(monos.size());
		auto perms = permutations_generator(rel_t::true_variables(variables));
		hasher<std::vector<int>> h(std::vector<int>(variables, 0), rel_t::max_exponent(variables, degree));
		for (const auto& i : monos) {
			auto hashed = h.hashvector(i.exponent);
			auto it = std::find(already_done.begin(), already_done.end(), hashed);
			if (it == already_done.end()) {
				auto o = orbit(i, perms);
				orbits.push_back(o);
				for (const auto& j : o)
					already_done.push_back(h.hashvector(j.exponent));
			}
		}
		return orbits;
	}

	///Class of elementary symmetric polynomials; no relations
	template<typename scalar_t>
	class elementary_symmetric {
	public:
		///Given vector exponent computes e_1^{exponent_1}*...*e_n^{exponent_n}
		polynomial<scalar_t, norelations> compute_product(const std::vector<int>& exponent) {
			polynomial<scalar_t, norelations> product = polynomial<scalar_t, norelations>::constant(1, elementary_symmetric_polynomials.size());
			for (int i = 0; i < exponent.size(); i++)
				if (exponent[i] != 0)
					product = multiply<scalar_t, norelations>(product, power(elementary_symmetric_polynomials[i], exponent[i], 0));
			return product;
		}

	protected:
		int n; ///<Number of Variables
		std::vector<polynomial<scalar_t, norelations>> elementary_symmetric_polynomials; ///<Vector containing the elementary symmetric polynomials e_1,...,e_n
		std::vector<std::string> elementary_symmetric_names; ///<Vector containing the names "e_i" of the elementary symmetric polynomials e_1,...,e_n

		///Constructs vector of elementary symmetric polynomials in terms of the variables x_1,...,x_n 
		elementary_symmetric(int n): n(n) {
			elementary_symmetric_polynomials.reserve(n);
			elementary_symmetric_names.reserve(n);
			for (int i = 1; i <= n; i++) {
				elementary_symmetric_polynomials.push_back(get_elementary_symmetric(i));
				elementary_symmetric_names.push_back("e_" + std::to_string(i));
			}
		}

	private:
		polynomial<scalar_t, norelations> get_elementary_symmetric(int i) {
			std::vector< monomial<scalar_t, norelations>> v;
			v.reserve(binomial(n, i));
			std::vector<int> mono(n);
			auto c = combinations_generator(n, i);
			for (const auto& comb : c) {
				for (const auto j : comb)
					mono[j] = 1;
				v.push_back(monomial<scalar_t, norelations>(1, mono));
				for (const auto j : comb)
					mono[j] = 0;
			}
			return polynomial<scalar_t, norelations>(v);
		}


	};

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	///Class for writing a symmetric polynomial as a polynomial on the elementary symmetric polynomials (a "decomposition")
	//
	///We denote the elementary symmetric polynomials by e_i
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	template<typename scalar_t>
	struct decomposition_elementary_symmetric : public elementary_symmetric<scalar_t> {
		polynomial<scalar_t, norelations> decomposition; ///< The polynomial on the elementary symmetric polynomials.

		///Initialize with number of variables
		decomposition_elementary_symmetric(int n) : elementary_symmetric<scalar_t>(n) {}

		///Initialize and decompose into polynomial on elementary symmetric polynomials
		decomposition_elementary_symmetric(const polynomial<scalar_t, norelations>& a) : elementary_symmetric<scalar_t>(a.monos[0].exponent.size()) { decompose(a); }


		///Decompose into polynomial on elementary symmetric polynomials
		void decompose(const polynomial<scalar_t, norelations>& a) {
			decomposition.monos.clear();
			decompose_recursive(a);
		}

		///Print a polynomial on elementary symmetric polynomials
		std::string print() const {
			return decomposition.print(this->elementary_symmetric_names);
		}

	private:
		void find_exponent(const monomial<scalar_t, norelations>& term, std::vector<int>& exponent) {
			for (int i = 0; i < term.exponent.size(); i++) {
				if (i == term.exponent.size() - 1)
					exponent.push_back(term.exponent[i]);
				else
					exponent.push_back(term.exponent[i] - term.exponent[i + 1]);
			}
		}
		void decompose_recursive(const polynomial<scalar_t, norelations>& a) {
			std::vector<int> pwr_here;
			pwr_here.reserve(a.monos[0].exponent.size());
			auto max = a.highest_term();
			find_exponent(max, pwr_here);
			auto product = this->compute_product(pwr_here);
			auto coeff = max.coeff / product.highest_term().coeff;
			decomposition.monos.push_back(monomial<scalar_t, norelations>(coeff, pwr_here));
			auto coeff_product = coeff * product;
			if (a != coeff_product)
				decompose_recursive(a - coeff_product);
		}


	};

	///Prints decomposition into elementary symmetric polynomials
	template<typename scalar_t>
	std::ostream& operator<<(std::ostream& os, const decomposition_elementary_symmetric<scalar_t>& a) {
		os << a.print();
		return os;
	}
}
