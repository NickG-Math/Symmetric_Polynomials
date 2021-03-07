#pragma once
#include <functional>
#include "General.h"

///@file
///@brief Contains classes for generating permutations, combinations and factory for such classes

namespace Symmetric_Polynomials {

	///////////////////////////////////////////////////////////////////////////////////
	///Prototype for classes providing iterators that generate elements such as interpolating vectors, permutations, combinations and so on
	//
	///Inherit from this class and define a method update() to get a const iterator (you will also need begin() and end() methods constructing such iterators; end() should always be defined by calling the factory end() ).
	///The template specialization is used for compile-time polymorphism (CRTP) so set it to be the child class.
	///The template gen_t is the type of the generated element.
	//////////////////////////////////////////////////////////////////////////////////
	template<typename specialization, typename gen_t>
	class Factory_Generator {
	public:
		///Returns the generated element
		const gen_t& operator*() const {
			return generated;
		}
		///Equality of iterators
		bool operator!=(const Factory_Generator& other) const {
			if (other.completed) //quick check if other=end
				return !completed;
			else
				return generated == other.generated;
		}
		///Generates next element
		specialization& operator++() {
			static_cast<specialization*>(this)->update();
			return *static_cast<specialization*>(this);
		}
		///Terminal iterator
		static specialization end() {
			specialization end;
			end.completed = 1;
			return end;
		}
	protected:
		///Default constructor
		Factory_Generator() {};
		bool completed; ///<1 if the iterator is end() i.e. if all elements have been generated
		gen_t generated; ///<The currently generated element
	};


	////////////////////////////////////////////////////////////////////////////////////////
	///Generates all permutations on a number of letters
	//
	//Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a permutations_generator object. Then i will be a permutation
	//////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	class Permutation_Generator {
		const T n;
	public:
		///Returns the total number of permutations i.e. the factorial
		size_t size() const {
			size_t factorial = 1;
			for (int i = 2; i <= n; i++) {
				factorial *= i;
			}
			return factorial;
		}
		///Constructor using number of letters
		Permutation_Generator(T n) : n(n) {}
		///Constant iterator that is used in a ranged for loop to generate the permutations. Non constant version is illegal
		class const_iterator : public Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>> {
			void update() {
				if (!std::next_permutation(this->generated.begin(), this->generated.end()))
					this->completed = 1;
			}
			const_iterator() {}
			friend class Permutation_Generator;
			friend class Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>;
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
			return Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>::end();
		}

	};



	////////////////////////////////////////////////////////////////////////////////////////
	///Generates all combinations on a number of letters making a number of choices
	//
	//Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a combinations_generator object. Then i will be a combination
	//////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	class Combination_Generator {
		const T total, choices;
	public:
		///Total number of combinations
		auto size() const {
			//safe implementation
			auto lowchoices = (choices > total - choices) ? total - choices : choices;
			auto binom = total;
			for (int i = 1; i <= lowchoices - 1; i++) { // $n/1 * (n-1)/2 * \cdots (n-k+1)/k$
				binom *= (total - i);
				binom /= i + 1;
			}
			return binom;
		}
		///Constructor using the total number of elements and the amount of choices we make
		Combination_Generator(T total, T choices) : total(total), choices(choices) {
			if (choices> total)
				throw("You can't choose more elements than those existing!");		
		}

		///Constant iterator that is used in a ranged for loop to generate the combinations. Non constant version is illegal
		class const_iterator : public Factory_Generator<Combination_Generator::const_iterator, std::vector<T>> {
			void update() {
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
			const_iterator() {}
			friend class Combination_Generator;
			friend class Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>;
			T total, choices;
		};
		///Initial generator
		const_iterator begin() const {
			const_iterator it;
			it.generated.resize(choices);
			std::iota(it.generated.begin(), it.generated.end(), 0);
			it.total = total;
			it.choices = choices;
			it.completed = 0;
			return it;
		}
		///Terminal generator.
		const_iterator end() const {
			return Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>::end();
		}
	};


	///Returns vector of all permutations on n letters
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

	/*
	namespace {//other unneeded generators

	////////////////////////////////////
	///The policy our monomials should satisfy. Policy is a char and takes values 0,1,-1. 1 if satisfied, -1 if not satisfied and we can guarantee that all larger vectors also won't satisfy it, and 0 otherwise.
	//
	//The -1 value speeds up computations comparad to 0 but is optional.
	////////////////////////////////////
		namespace policy {

			///The always true policy
			template<typename T>
			static char always_true(const T&) {
				return 1;
			}

			///The policy made by a given comparactor function
			template<typename T>
			static std::function<char(const std::vector<char>&)> custom(const T& comparator) {
				return std::bind(&T::acceptable, comparator, std::placeholders::_1); //bind the function acceptable to this comparator instance. The placeholder is for the one vector argument
			}
		};

		////////////////////////////////////////////////////////////////////////////////////////
		///Generates all vectors interpolating between a given min and max vector (with the same length), and satisfying a given policy
		//
		//Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a vector_interpolate_generator object. Then i will be an interpolating vector satisfying the policy
		//////////////////////////////////////////////////////////////////////////////////
		class vector_interpolate_generator {
		public:
			///Constructor using min, max and policy
			vector_interpolate_generator(const std::vector<char>& min, const std::vector<char>& max, const std::function<char(const std::vector<char>&)>& policy = policy::always_true<std::vector<char>>)
				: min(min), max(max), policy(policy) {}
			///Returns upper bound on the amount of generated elements
			size_t size() {
				if (min.empty())
					return 0;
				size_t total = max[0] - min[0] + 1;
				for (int i = 1; i < max.size(); i++) {
					total *= max[i] - min[i] + 1;
				}
				return total;
			}

			///Constant iterator that is used in a ranged for loop to generate the interpolating vectors. Non constant version is illegal
			class const_iterator : public Factory_Generator<vector_interpolate_generator::const_iterator, std::vector<char>> {
				void get_first_satisfying_policy() {
					policy_status = (*policy)(generated);
					while (policy_status != 1 && !completed) //update until its satisfied
						update();
				}
				void update() {
					do {
						int i = 0;
						bool current_policy = 1;
						while (i < length && ((generated[i] == *(max_ptr + i) || (policy_status == -1 && current_policy)))) {
							current_policy = 0; //the overshoot may only happen once per loop, as the policy is only updated outside the loop
							generated[i] = *(min_ptr + i);
							i++;
						}
						if (i == length) {
							completed = 1;
							return;
						}
						generated[i]++;
						policy_status = (*policy)(generated);
					} while (policy_status != 1 && !completed);
				}
				const_iterator() {};
				char policy_status;
				int length;
				const char* min_ptr;
				const char* max_ptr;
				const std::function<char(const std::vector<char>&)>* policy;
				friend class vector_interpolate_generator;
				friend class Factory_Generator<vector_interpolate_generator::const_iterator, std::vector<char>>;

			};
			///Initial generator
			const_iterator begin() const {
				const_iterator it;
				it.length = min.size();
				it.generated = min;
				it.max_ptr = max.data();
				it.min_ptr = min.data();
				it.policy = &policy;
				it.get_first_satisfying_policy(); //in case the min does not satisfy the policy
				it.completed = 0;
				return it;
			}
			///Terminal generator.
			const_iterator end() const {
				return Factory_Generator<vector_interpolate_generator::const_iterator, std::vector<char>>::end();
			}
		private:
			const std::vector<char> min, max;
			const std::function<char(const std::vector<char>&)> policy;
		};



		///Returns all vectors interpolating between min, max and satisfying a given policy
		std::vector<std::vector<char>> vector_interpolate(const std::vector<char>& min, const std::vector<char>& max, const std::function<char(const std::vector<char>&)>& policy = policy::always_true<std::vector<char>>) {
			std::vector<std::vector<char>> vectors;
			vector_interpolate_generator gen(min, max, policy);
			vectors.reserve(gen.size());
			for (const auto& i : gen)
				vectors.push_back(i);
			return vectors;
		}


		///Applies permutation perm on the target from starting point begin
		template<typename Iterator>
		void apply_permutation(std::vector<char>& target, Iterator begin, const std::vector<char>& perm) {
			for (const auto i : perm)
				target.push_back(*(begin + i));
		}


		///Breaks vector v into given number of pieces and applies given permutation on each one; then joins the pieces and returns the vector
		std::vector<char> apply_permutation_pieces(const std::vector<char>& v, int pieces, const std::vector<char>& perm) {
			std::vector<char> target;
			target.reserve(v.size());
			for (int i = 0; i < pieces; i++)
				apply_permutation(target, v.begin() + i * perm.size(), perm);
			return target;
		}

		///Applies permutation perm on the vector v
		std::vector<char> apply_permutation(const std::vector<char>& v, const std::vector<char>& perm) {
			return apply_permutation_pieces(v, 1, perm);
		}
	}
	*/
}
