#pragma once
#include "../Generators.hpp"
///@file
///@brief Contains classes for generating permutations, combinations and factory for such classes

namespace Symmetric_Polynomials {

	template<typename specialization, typename gen_t>
	const gen_t& Factory_Generator<specialization, gen_t>::operator*() const {
		return generated;
	}

	template<typename specialization, typename gen_t>
	bool Factory_Generator<specialization, gen_t>::operator!=(const Factory_Generator& other) const {
		if (other.completed) //quick check if other=end
			return !completed;
		else
			return generated == other.generated;
	}

	template<typename specialization, typename gen_t>
	specialization& Factory_Generator<specialization, gen_t>::operator++() {
		static_cast<specialization*>(this)->update();
		return *static_cast<specialization*>(this);
	}

	template<typename specialization, typename gen_t>
	specialization Factory_Generator<specialization, gen_t>::end() {
		specialization end;
		end.completed = 1;
		return end;
	}


	template<typename T>
	Permutation_Generator<T>::Permutation_Generator(T n) : n(n) {}

	template<typename T>
	size_t Permutation_Generator<T>::size() const {
		size_t factorial = 1;
		for (int i = 2; i <= n; i++) {
			factorial *= i;
		}
		return factorial;
	}

	template<typename T>
	void Permutation_Generator<T>::const_iterator::update() {
		if (!std::next_permutation(this->generated.begin(), this->generated.end()))
			this->completed = 1;
	}

	template<typename T>
	typename Permutation_Generator<T>::const_iterator Permutation_Generator<T>::begin() const {
		const_iterator it;
		it.generated.resize(n);
		std::iota(it.generated.begin(), it.generated.end(), 0);
		it.completed = 0;
		return it;
	}

	template<typename T>
	typename Permutation_Generator<T>::const_iterator Permutation_Generator<T>::end() const {
		return Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>::end();
	}


	template<typename T>
	Combination_Generator<T>::Combination_Generator(T total, T choices) : total(total), choices(choices) {
		if (choices > total) {
			std::cerr << "You can't choose more elements than those existing!";
			abort();
		}
	}

	template<typename T>
	auto Combination_Generator<T>::size() const {
		//safe implementation
		auto lowchoices = (choices > total - choices) ? total - choices : choices;
		auto binom = total;
		for (int i = 1; i <= lowchoices - 1; i++) { // $n/1 * (n-1)/2 * \cdots (n-k+1)/k$
			binom *= (total - i);
			binom /= i + 1;
		}
		return binom;
	}


	template<typename T>
	void Combination_Generator<T>::const_iterator::update() {
		for (int i = (int)choices - 1; i >= 0; i--) { //must be signed!
			if (this->generated[i] < total - choices + i) {
				this->generated[i]++;
				for (auto j = i + 1; j < (int)choices; j++)
					this->generated[j] = this->generated[j - 1] + 1;
				return;
			}
		}
		this->completed = 1;
		return;
	}

	template<typename T>
	typename Combination_Generator<T>::const_iterator Combination_Generator<T>::begin() const {
		const_iterator it;
		it.generated.resize(choices);
		std::iota(it.generated.begin(), it.generated.end(), 0);
		it.total = total;
		it.choices = choices;
		it.completed = 0;
		return it;
	}

	template<typename T>
	typename Combination_Generator<T>::const_iterator Combination_Generator<T>::end() const {
		return Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>::end();
	}

	template<typename T>
	std::vector<std::vector<T>> all_permutations(T n) {
		std::vector<std::vector<T>> v;
		Permutation_Generator<T> p(n);
		v.reserve(p.size());
		for (const auto& i : p)
			v.push_back(i);
		return v;
	}

	///Returns vector of all combinations on n letters choosing m many  
	template<typename T>
	std::vector<std::vector<T>> all_combinations(T n, T m) {
		std::vector<std::vector<T>> v;
		Combination_Generator<T> p(n, m);
		v.reserve(p.size());
		for (const auto& i : p)
			v.push_back(i);
		return v;
	}
}
