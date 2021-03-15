#pragma once
#include <algorithm>
#include "General.hpp"

///@file
///@brief Contains classes for generating permutations, combinations and factory for such classes

namespace Symmetric_Polynomials {

	///////////////////////////////////////////////////////////////////////////////////
	///Prototype for classes providing iterators that generate elements such as interpolating vectors, permutations, combinations and so on
	//
	///Inherit from this class and define a method update() to get a const iterator.
	///You will also need begin() and end() methods constructing such iterators; end() should always be defined by calling the factory end().
	///Warning: Probably not thread-safe (depends on child's update() method).
	///@tparam specialization Used for compile-time polymorphism (CRTP): set it to be the child class.
	///@tparam gen_t The type of the generated element.
	//////////////////////////////////////////////////////////////////////////////////
	template<typename specialization, typename gen_t>
	class Factory_Generator {
	public:
		///Returns the generated element
		const gen_t& operator*() const;
		///Equality of iterators
		bool operator!=(const Factory_Generator& other) const;
		///Generates next element
		specialization& operator++();
		///Terminal iterator
		static specialization end();
	protected:
		///Default constructor
		Factory_Generator()=default;
		bool completed; ///<1 if the iterator is end() i.e. if all elements have been generated
		gen_t generated; ///<The currently generated element
	};


	////////////////////////////////////////////////////////////////////////////////////////
	///Generates all permutations on a number of letters
	//
	///Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a permutations_generator object. Then i will be a permutation
	///Warning: Not thread safe!
	///@tparam T The data type of our permutations eg std::vector<int>
	//////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	class Permutation_Generator {
		const T n;
	public:
		
		///Computes total number of permutations
		//
		///@return \f$n!\f$ where \f$n\f$ is the number of letters
		size_t size() const;

		///Sets up the generator
		//
		///@param n the total number of letters
		Permutation_Generator(T n);

		///Constant iterator that is used in a ranged for loop to generate the permutations. Non constant version is illegal
		class const_iterator : public Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>> {
			void update();
			const_iterator() = default;
			friend class Permutation_Generator;
			friend class Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>;
		};
		///Initial generator
		const_iterator begin() const;
		///Terminal generator.
		const_iterator end() const;
	};



	////////////////////////////////////////////////////////////////////////////////////////
	///Generates all combinations on a number of letters making a number of choices
	//
	///Use with a ranged for loop: Ex for (const auto& i:v){ ... } where v is a combinations_generator object. Then i will be a combination
	///Warning: Not thread safe!
	///@tparam T The data type of our combinations eg std::vector<int>
	//////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	class Combination_Generator {
		const T total, choices;
	public:
		///Computes total number of combinations
		//
		///@return \f${n}\choose{k}\f$ where \f$n\f$=total and \f$k\f$=choices
		auto size() const;

		///Sets up the generator
		//
		///@param total the total number of letters
		///@param choices the total number of choices
		Combination_Generator(T total, T choices);

		///Constant iterator that is used in a ranged for loop to generate the combinations. Non constant version is illegal
		class const_iterator : public Factory_Generator<Combination_Generator::const_iterator, std::vector<T>> {
			void update();
			const_iterator() = default;
			friend class Combination_Generator;
			friend class Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>;
			T total, choices;
		};
		///Initial generator
		const_iterator begin() const;
		///Terminal generator.
		const_iterator end() const;
	};


	///Returns vector of all permutations on n letters
	//
	///@tparam T The value type of the permutations
	///@param n The number of letters
	///@return Vector of vectors each of which is \f$[a_1,...,a_n]\f$ with \f$0\le a_i\le n\f$ all distinct
	template<typename T>
	std::vector<std::vector<T>> all_permutations(T n);

	///Returns vector of all combinations on n letters choosing m many  
	//
	///@tparam T The value type of the combinations
	///@param n The total number of letters
	///@param m The number of choices 
	///@return Vector of vectors each of which is \f$[a_1,...,a_n]\f$ with \f$0\le a_i\le n\f$ all distinct
	template<typename T>
	std::vector<std::vector<T>> all_combinations(T n, T m);
}
#include "impl/Generators.ipp"