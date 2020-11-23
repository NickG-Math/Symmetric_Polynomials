#pragma once
#include "Relations.h"

///@file
///@brief Contains generators of all monomials with desired degree and satisfying a desired policy.

namespace Symmetric_Polynomials {

	///////////////////////////////////////////////////////////////////////////////////
	///Prototype for an iterator that generates elements such as interpolating vectors, permutations, combinations and so on
	//
	//Inherit from this class and define a method update() to get a const iterator (you will also need begin() and end() methods constructing such iterators; end() should always be defined by calling returning terminal() ).
	//The template specialization is used for compile-time polymorphism (CRTP). Set it to be the child class and make sure the child class befriends generating_const_iterator.
	//The template gen_t is the type of the generated element.
	template<typename specialization, typename gen_t>
	class generating_const_iterator {
	public:
		///Returns the generated element
		const gen_t& operator*() const {
			return generated;
		}
		///Equality of iterators
		bool operator!=(const generating_const_iterator& other) const {
			if (other.completed) //quick check
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
		static specialization terminal() {
			specialization end;
			end.completed = 1;
			return end;
		}
	protected:
		///Default constructor
		generating_const_iterator() {};
		bool completed; ///<1 if the iterator is end() i.e. if all elements have been generated
		gen_t generated; ///<The currently generated element
	};



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
		static std::function<char(const std::vector<int>&)> custom(const T& comparator) {
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
		vector_interpolate_generator(const std::vector<int>& min, const std::vector<int>& max, const std::function<char(const std::vector<int>&)>& policy = policy::always_true<std::vector<int>>)
			: min(min), max(max), policy(policy) {}
		///Returns upper bound on the amount of generated elements
		long size() {
			if (min.empty())
				return 0;
			long total = max[0] - min[0] + 1;
			for (int i = 1; i < max.size(); i++) {
				total *= max[i] - min[i] + 1;
			}
			return total;
		}
			   
		///Constant iterator that is used in a ranged for loop to generate the interpolating vectors. Non constant version is illegal
		class const_iterator : public generating_const_iterator<vector_interpolate_generator::const_iterator, std::vector<int>> {
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
			const int* min_ptr;
			const int* max_ptr;
			const std::function<char(const std::vector<int>&)>* policy;
			friend class vector_interpolate_generator;
			friend class generating_const_iterator<vector_interpolate_generator::const_iterator, std::vector<int>>;

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
			it.completed=0;
			return it;
		}
		///Terminal generator.
		const_iterator end() const {
			return generating_const_iterator<vector_interpolate_generator::const_iterator, std::vector<int>>::terminal();
		}
	private:
		const std::vector<int> min, max;
		const std::function<char(const std::vector<int>&)> policy;
	};



	///Returns all vectors interpolating between min, max and satisfying a given policy
	std::vector<std::vector<int>> vector_interpolate(const std::vector<int>& min, const std::vector<int>& max, const std::function<char(const std::vector<int>&)>& policy = policy::always_true<std::vector<int>>) {
		std::vector<std::vector<int>> vectors;
		vector_interpolate_generator gen(min, max, policy);
		vectors.reserve(gen.size());
		for (const auto& i : gen)
			vectors.push_back(i);
		return vectors;
	}

	///The comparator that checks if the degree of a given exponent vector is equal to a desired degree
	template<typename rel_t>
	struct degree_comparator {
		///Returns 1 if the degree agrees with the desired degree, -1 if it's higher than the desired degree and 0 otherwise.
		char acceptable(const std::vector<int>& a) {
			auto d = rel_t::compute_degree(a);
			if (d > desired_degree)
				return -1;
			else if (d == desired_degree)
				return 1;
			else
				return 0;
		}
		///Constructor given the desired degree
		degree_comparator(int desired_degree) : desired_degree(desired_degree) {}

	private:
		const int desired_degree;
	};

	//////////////////////////////////////////////////////////////////
	///The comparator that checks  if the degree of a given exponent vector is equal to a desired degree, and further disqualifies a exponent vector if it contains a relation.
	//
	///Eg a relation can be \f$a_1^2a_2=\cdots\f$. In that case, an exponent vector [2,1,0,...,] is disqualified, and so are all [v_1,v_2,v_3,...] with v_1>=2 and v_2>=1 since these also contain the relation
	///The relations are allowed to be very complicated (think equivariant Chern classes) hence why they are not encoded in a rel_t
	/////////////////////////////////////////////////////////////////////////////
	class degree_comparator_relations {
		const int desired_degree;
		const std::vector<int> dimensions;
		const std::vector<std::vector<char>> relations;

		bool has_relation(const std::vector<int>& a) {
			for (const auto& rel : relations)
				if (rel <= a)
					return 1;
			return 0;
		}

	public:
		///Constructor given the desired degree, unacceptable combinations and dimensions
		degree_comparator_relations(int desired_degree, const std::vector<std::vector<char>>& relations, const std::vector<int>& dimensions)
			: desired_degree(desired_degree), dimensions(dimensions) ,relations(relations) {}


		///Returns 1 if the degree agrees with the desired degree and has no relations, -1 if it's higher than the desired degree and 0 otherwise.
		char acceptable(const std::vector<int>& a) {
			auto d = general_compute_degree(a, dimensions);
			if (d > desired_degree)
				return -1;
			else if (d < desired_degree)
				return 0;
			return !has_relation(a);
		}
	};


	///Returns the monomial basis in given number of variables and degree
	template<typename scalar_t, typename rel_t>
	std::vector<monomial<scalar_t, rel_t>> monomial_basis(int variables, int degree) {
		std::vector<monomial<scalar_t, rel_t>> monos;
		vector_interpolate_generator gen(std::vector<int>(variables, 0), rel_t::max_exponent(variables, degree), policy::custom<degree_comparator<rel_t>>(degree_comparator<rel_t>(degree)));
		monos.reserve(gen.size());
		for (const auto& i : gen)
			monos.push_back(monomial<scalar_t, rel_t>(1, i));
		return monos;
	}
}
