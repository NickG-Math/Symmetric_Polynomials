#pragma once
#include "Polynomials.h"

///@file
///@brief Contains methods for writing symmetric polynomials (with no relations) in terms of the elementary symmetric ones

namespace Symmetric_Polynomials {

	///Variables \f$e_1,...,e_n\f$ denoting the elementary symmetric polynomials \f$e_i=\sigma_i\f$, of degrees \f$|e_i|=i\f$
	template<typename T>
	struct Elementary_Symmetric_Variables : public Standard_Variables<T> {
		using Standard_Variables<T>::Standard_Variables;
		///The degree of a monomial on the \f$e_i\f$, \f$|e_i|=i\f$
		int degree() const {
			int deg = 0;
			for (int i = 0; i < size(); i++)
				deg += (*this)[i] * (i + 1);
			return deg;
		}
		///The names of the \f$e_i\f$
		static std::string name(int i) {
			return "e_" + std::to_string(i + 1);
		}
	};



	///Class for symmetric polynomials with no relations, allowing transformation from \f$x_i\f$ variables to \f$e_i\f$ variables and vice-versa.
	template<typename scalar_t, typename exponent_t>
	class Symmetric_Basis {
		typedef Elementary_Symmetric_Variables<typename exponent_t::value_type> elementary_symmetric_t;
	public:

		///Transform a polynomial on the elementary symmetric polynomials e_i into a symmetric polynomial on the original variables x_i
		Polynomial<scalar_t, exponent_t> operator()(const Polynomial<scalar_t, elementary_symmetric_t>& poly_on_ei) const {
			Polynomial<scalar_t, exponent_t> p(n);
			for (auto it = poly_on_ei.begin(); it != poly_on_ei.end(); ++it) {
				auto prod = compute_product(it.exponent());
				prod *= it.coeff();
				p += prod;
			}
			return p;
		}

		///Transform a symmetric polynomial on the original variables x_i into one on the elementary symmetric polynomials e_i
		Polynomial<scalar_t, elementary_symmetric_t> operator()(Polynomial<scalar_t, exponent_t> a) const  {
			Polynomial<scalar_t, elementary_symmetric_t> decomposition(n);
			elementary_symmetric_t pwr_here;
			pwr_here.reserve(n);
			while (true) {
				pwr_here.clear();
				auto max = a.highest_term();
				find_exponent(max.exponent, pwr_here);
				auto product = this->compute_product(pwr_here);
				auto coeff = max.coeff / product.highest_term().coeff;
				decomposition.insert(pwr_here, coeff);
				product *= coeff;
				if (a == product)
					return decomposition;
				else
					a -= product;
			}
		}

		///Constructor given number of variables
		Symmetric_Basis(int n) : n(n) {
			elementary_symmetric_polynomials.reserve(n);
			for (int i = 1; i <= n; i++)
				elementary_symmetric_polynomials.push_back(get_elementary_symmetric(i));
		}
	private:

		const int n; ///<Number of Variables
		std::vector<Polynomial<scalar_t, exponent_t>> elementary_symmetric_polynomials; ///<Vector containing the elementary symmetric polynomials e_1,...,e_n


		Polynomial<scalar_t, exponent_t> get_elementary_symmetric(int i) const {
			Polynomial<scalar_t, exponent_t> poly(n);
			exponent_t mono(n);
			auto c = Combination_Generator(n, i);
			for (const auto& comb : c) {
				for (const auto j : comb)
					mono[j] = 1;
				poly.insert(mono, 1);
				for (const auto j : comb)
					mono[j] = 0;
			}
			return poly;
		}

		void find_exponent(const exponent_t& term, exponent_t& exponent) const {
			for (int i = 0; i < term.size(); i++) {
				if (i == term.size() - 1)
					exponent.push_back(term[i]);
				else
					exponent.push_back(term[i] - term[i + 1]);
			}
		}

		///Given vector exponent computes e_1^{exponent_1}*...*e_n^{exponent_n}
		Polynomial<scalar_t, exponent_t> compute_product(const elementary_symmetric_t& exponent) const {
			auto product = Polynomial<scalar_t, exponent_t>::constant(1, elementary_symmetric_polynomials.size());
			for (int i = 0; i < exponent.size(); i++)
				if (exponent[i] != 0)
					product *= (elementary_symmetric_polynomials[i] ^ exponent[i]);
			return product;
		}


	};

}
