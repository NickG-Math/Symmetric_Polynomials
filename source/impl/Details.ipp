#pragma once
#include <type_traits>
#include <utility>
#include "stddef.h"

///	@file
///	@brief Contains various implementation details

#if defined SYMMP_USE_OPEN_MP & defined _OPENMP
///	@brief		Macro that parallelizes certain loops via openMP if SYMMP_USE_OPEN_MP is defined, and does nothing otherwise
///	@details	Only used for the loops in \c print_half_idempotent_relations and \c pontryagin_via_chern \n
///				To configure the number of threads please set the environment variable OMP_NUM_THREADS on the console before running any executable
#define SYMMP_RUN_LOOP_IN_PARALLEL _Pragma("omp parallel for schedule(dynamic)")
#if defined _MSC_VER
#pragma message("symmp: openMP enabled!")
#else
#pragma message "symmp: openMP enabled!"
#endif
#else
///	@brief		Macro that parallelizes certain loops via openMP if SYMMP_USE_OPEN_MP is defined, and does nothing otherwise
///	@details	Only used for the loops in \c print_half_idempotent_relations and \c pontryagin_via_chern \n
///				To configure the number of threads please set the environment variable OMP_NUM_THREADS on the console before running any executable
#define SYMMP_RUN_LOOP_IN_PARALLEL
#if defined _MSC_VER
#pragma message("symmp: openMP not enabled! To enable please define SYMMP_USE_OPEN_MP before including any header of this library and use the relevant compiler option eg -fopenmp")
#else
#pragma message "symmp: openMP not enabled! To enable please define SYMMP_USE_OPEN_MP before including any header of this library and use the relevant compiler option eg -fopenmp"
#endif
#endif


namespace symmp
{

	///Contains various implementation details such as SFINAE
	namespace implementation_details {

		template <typename _exp>
		using pair_t = std::pair<typename _exp::deg_t, _exp>;

		///Given a pair, hash only the second parameter
		template <typename _exp>
		struct hash_only_exp
		{
			///Hash only second parameter
			auto operator()(const pair_t<_exp>& pair) const
			{
				return pair.second();
			}

			auto operator()(const _exp& a) const {
				return a();
			}
		};

		template <class _scl, class _exp, template<class...> class _cnt, bool _ord, class ... _arg>
		using container_wrapper = std::conditional_t<_ord, _cnt<pair_t<_exp>, _scl, _arg...>, _cnt<pair_t<_exp>, _scl, hash_only_exp<_exp>, _arg...>>;

		template <typename T>
		static constexpr std::false_type test_degree_existence(...);

		template <typename T>
		static constexpr decltype(std::declval<T>().degree(), std::true_type()) test_degree_existence(int);

		template <typename T>
		using has_degree_function = decltype(test_degree_existence<T>(0));

		template <typename T>
		static constexpr std::false_type test_name_existence(...);

		template <typename T>
		static constexpr decltype(std::declval<T>().name(std::declval<int>(), std::declval<int>()), std::true_type()) test_name_existence(int);

		template <typename T>
		using has_name_function = decltype(test_name_existence<T>(0));
	}
}