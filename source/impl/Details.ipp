#pragma once
#include <type_traits>

namespace symmp
{

	namespace implementation_details {

		template <typename _exp>
		using pair_t = std::pair<typename _exp::deg_t, _exp>;

		///Given a pair, hash only the second parameter
		template <typename _exp>
		struct hash_only_exponent
		{
			///Hash only second parameter
			auto operator()(const pair_t<_exp>& pair) const
			{
				return pair.second();
			}
		};

		template <typename _scl, typename _exp, template<typename...> typename _container, bool container_is_ordered, typename ... _Args>
		struct polynomial_data_type;

		template <typename _scl, typename _exp, template<typename...> typename _container, typename ... _Args>
		struct polynomial_data_type<_scl, _exp, _container, 1, _Args...>: _container<pair_t<_exp>, _scl, _Args...> {
			using _container<pair_t<_exp>, _scl, _Args...>::_container;
			friend class Polynomial;
		};

		template <typename _scl, typename _exp, template<typename...> typename _container, typename ... _Args>
		struct polynomial_data_type<_scl, _exp, _container, 0, _Args...>: _container<pair_t<_exp>, _scl, hash_only_exponent<_exp>, _Args...> {
			using _container<pair_t<_exp>, _scl, hash_only_exponent<_exp>, _Args...>::_container;
			friend class Polynomial;
		};


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

		template <typename T>
		static constexpr std::false_type test_reserve_existence(...);

		template <typename T>
		static constexpr decltype(std::declval<T>().reserve(std::declval<size_t>()), std::true_type()) test_reserve_existence(int);

		template <typename T>
		using has_reserve_function = decltype(test_reserve_existence<T>(0));

		template <typename T, typename S>
		std::pair<T, S> operator+(const std::pair<T, S>& a, const std::pair<T, S> b)
		{
			return std::pair<T, S>(a.first + b.first, a.second + b.second);
		}
	}
}