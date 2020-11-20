#pragma once
#include "Relations.h"

///@file
///@brief Contains generators of all monomials with desired degree and satisfying a desired policy.

namespace Symmetric_Polynomials{

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
///Generates all vectors between a given min and max vector, and satisfying a given policy
//
///Don't use with empty vectors with this!
//////////////////////////////////////////////////////////////////////////////////
class vector_interpolate_generator {
public:
	long total; ///<The total number of vectors

	///Constructor using min, max and policy
	vector_interpolate_generator(const std::vector<int>& min, const std::vector<int>& max, const std::function<char(const std::vector<int>&)>& policy = policy::always_true<std::vector<int>>)
		: min(min), max(max), policy(policy) {
		total = max[0] - min[0] + 1;
		for (int i = 1; i < max.size(); i++) {
			total *= max[i] - min[i] + 1;
		}
		initialize();
	}

	///Returns the generated vector
	std::vector<int> get_generated() const {
		return generated;
	}

	///False if all vectors have been generated
	inline bool keep_going() {
		return !completed;
	}

	///Initialize generator
	void initialize(){
		completed = 0;
		generated = min;
		policy_status = policy(generated);
	}

	///Checks if policy is satisfied
	bool policy_satisfied() {
		return (policy_status == 1); //policy_status is char
	}

	///Generates new vector
	void update() {
		int i = 0;
		bool current_policy=1;
		while (i < max.size() && ((generated[i] == max[i] || (policy_status == -1 && current_policy)))) {
			current_policy = 0; //the overshoot may only happen once per loop, as the policy is only updated outside the loop
			generated[i] = min[i];
			i++;
		}
		if (i == max.size()) {
			completed=1;
			return;
		}
		generated[i]++;
		policy_status = policy(generated);
	}
private:
	const std::vector<int> min, max;
	const std::function<char(const std::vector<int>&)> policy;
	bool completed;
	char policy_status;
	std::vector<int> generated;
};

///Returns all vectors between min, max and satisfying a given policy
std::vector<std::vector<int>> vector_interpolate(const std::vector<int>& min, const std::vector<int>& max, const std::function<char(const std::vector<int>&)>& policy = policy::always_true<std::vector<int>>) {
	std::vector<std::vector<int>> vectors;
	if (min.empty()) {
		if (policy(min))
			vectors.push_back(min);
		return vectors;
	}
	vector_interpolate_generator gen(min, max, policy);
	vectors.reserve(gen.total);
	while (gen.keep_going()) {
		if (gen.policy_satisfied())
			vectors.push_back(gen.get_generated());
		gen.update();
	}
	return vectors;
}

///The comparator that checks if the degree of a given powers vector is equal to a desired degree
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
///The comparator that checks  if the degree of a given powers vector is equal to a desired degree, and further disqualifies a powers vector if it contains a relation.
//
///Eg a relation can be \f$a_1^2a_2=\cdots\f$. In that case, a power vector [2,1,0,...,] is disqualified, and so are all [v_1,v_2,v_3,...] with v_1>=2 and v_2>=1 since these also contain the relation
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
		: desired_degree(desired_degree), relations(relations), dimensions(dimensions) {}


	///Returns 1 if the degree agrees with the desired degree and has no relations, -1 if it's higher than the desired degree and 0 otherwise.
	char acceptable(const std::vector<int>& a) {
		auto d = general_compute_degree(a,dimensions);
		if (d > desired_degree)
			return -1;
		else if (d < desired_degree)
			return 0;
		return !has_relation(a);
	}
};

///Generates a basis of monomials in given number of variables and degree
template<typename scalar_t, typename rel_t>
struct monomial_basis_generator : public vector_interpolate_generator {
	///Constructor given number of variables and degree
	monomial_basis_generator(int variables, int degree)
		: vector_interpolate_generator(std::vector<int>(variables, 0), rel_t::max_powers(variables, degree), policy::custom<degree_comparator<rel_t>>(degree_comparator<rel_t>(degree))) {}
	///Returns the generated monomial
	monomial<scalar_t, rel_t> get_generated() const {
		return monomial<scalar_t, rel_t>(1, vector_interpolate_generator::get_generated());
	}
};

///Returns the monomial basis in given number of variables and degree
template<typename scalar_t, typename rel_t>
std::vector<monomial<scalar_t, rel_t>> monomial_basis(int variables, int degree) {
	std::vector<monomial<scalar_t, rel_t>> monos;
	monomial_basis_generator<scalar_t, rel_t>gen(variables, degree);
	monos.reserve(gen.total);
	while (gen.keep_going()) {
		if (gen.policy_satisfied())
			monos.push_back(gen.get_generated());
		gen.update();
	}
	return monos;
}
}
