#pragma once
#include <algorithm>
#include "General.hpp"

///@file
///@brief Contains classes for generating permutations, combinations and a factory for such classes

namespace Symmetric_Polynomials
{

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief			Prototype for coroutine-like iterators that generate elements such as interpolating vectors, permutations, combinations...
	///	@details		Inherit from this class and define a method update() to get a const iterator. \n
	///					You will also need begin() and end() methods constructing such iterators; end() should always be defined by calling the factory end().\n
	///					Example implementations: \c Combination_Generator and \c Permutation_Generator
	///	@attention		Probably not thread-safe (depends on child's update() method).
	///	@todo			Implement as coroutine (C++20)
	///	@tparam	spec_t	Used for compile-time polymorphism (CRTP): set it to be the child class.
	///	@tparam	gen_t	The type of the generated element.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename spec_t, typename gen_t>
	class Factory_Generator
	{
	public:
		const gen_t &operator*() const;						   ///<Returns the generated element
		bool operator!=(const Factory_Generator &other) const; ///<Inequality of iterators (used to detect if generation has been completed)
		spec_t &operator++();								   ///<Generates next element
		static spec_t end();								   ///<Terminal iterator
	protected:
		bool completed;	 ///<1 if the iterator is end() i.e. if all elements have been generated
		gen_t generated; ///<The currently generated element
	};

	////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Generates all permutations on a number of letters
	///	@details	Use with a ranged for loop: \code for (const auto& i:v) {...} \endcode where \c v is a \c Permutation_Generator object. Then \c i will be a permutation
	///	@warning	Not thread safe!
	///	@todo		Implement via coroutine (C++20)
	///	@tparam T	The data type of our permutations eg \c std::vector<int>
	//////////////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Permutation_Generator
	{
		const T n;

	public:
		/// @brief	Computes total number of permutations
		/// @return \f$n!\f$ where \f$n\f$ is the number of letters
		size_t size() const;

		///	@brief		Constructor sets up the generator
		///	@param n	The total number of letters
		Permutation_Generator(T n);

		///	@brief		Constant iterator that is used in a ranged for loop to generate the permutations.
		///	@warning	Non constant version is illegal
		class const_iterator : public Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>
		{
			void update();
			friend class Permutation_Generator;													   ///<Befriending outer class
			friend class Factory_Generator<Permutation_Generator::const_iterator, std::vector<T>>; ///<Befriending parent
		};

		/// @brief	Begin iterator
		/// @return An iterator to the first generated element
		const_iterator begin() const;

		/// @brief	End iterator
		/// @return An iterator to the end of the generator (equality with this indicates that the generator has completed)
		const_iterator end() const;
	};

	////////////////////////////////////////////////////////////////////////////////////////
	///	@brief		Generates all combinations on a number of letters making a number of choices
	///	@details	Use with a ranged for loop: \code for (const auto& i:v) {...} \endcode where \c v is a \c Combination_Generator object. Then \c i will be a combination
	///	@warning	Not thread safe!
	///	@todo		Implement via coroutine (C++20)
	///	@tparam  T	The data type of our combinations eg \c std::vector<int>
	//////////////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Combination_Generator
	{
		const T total, choices;

	public:
		///	@brief	Computes total number of combinations
		///	@return \f${n}\choose{k}\f$ where \f$n\f$=total and \f$k\f$=choices
		auto size() const;

		///	@brief			Sets up the generator
		///	@param total	The number of letters
		///	@param choices	The number of choices
		Combination_Generator(T total, T choices);

		///	@brief		Constant iterator that is used in a ranged for loop to generate the combinations.
		///	@warning	Non \c const version is illegal
		class const_iterator : public Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>
		{
			void update();
			friend class Combination_Generator;													   ///<Befriending outer class
			friend class Factory_Generator<Combination_Generator::const_iterator, std::vector<T>>; ///<Befriending parent
			T total, choices;
		};

		/// @brief	Begin iterator
		/// @return An iterator to the first generated element
		const_iterator begin() const;

		/// @brief	End iterator
		/// @return An iterator to the end of the generator (equality with this indicates that the generator has completed)
		const_iterator end() const;
	};

	///@brief		Returns vector of all permutations on \p n letters
	///@tparam	T	The value type of the permutations
	///@param	n	The number of letters
	///@return		Vector of vectors each of which is \f$[a_1,...,a_n]\f$ with \f$0\le a_i\le n\f$ all distinct
	template <typename T>
	std::vector<std::vector<T>> all_permutations(T n);

	///@brief		Returns vector of all combinations on \p n letters choosing \p m many
	///@tparam	T	The value type of the combinations
	///@param	n	The number of letters
	///@param	m	The number of choices
	///@return		Vector of vectors each of which is \f$[a_1,...,a_n]\f$ with \f$0\le a_i\le n\f$ all distinct
	template <typename T>
	std::vector<std::vector<T>> all_combinations(T n, T m);
}
#include "impl/Generators.ipp"