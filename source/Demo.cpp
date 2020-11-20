#include "Half_Idempotent.h"

///@file
///@brief A demonstration file that can be compiled


int main() {
	using namespace Symmetric_Polynomials;
	monomial<rational, norelations> m(rational(2, 3), { 0,1,4 });
	std::cout << m << "\n";
	auto u = m.print();
	auto v = m.print({ "x","y", "z" });
	polynomial<int, halfidempotent> p({ monomial<int,halfidempotent>(2,{1,2,0,1}), monomial<int,halfidempotent>(3,{3,0,1,1}) });
	decomposition_elementary_symmetric<int> dec(polynomial<int, norelations>({ monomial<int, norelations>(2, { 1,2 }), monomial<int, norelations>(2, {2,1 }) }));
	std::cout<< dec << "\n";
	decomposition_half_idempotent<rational> dec2(polynomial<rational, halfidempotent>({ monomial<rational, halfidempotent>(2, { 1,0,1,0 }), monomial<rational, halfidempotent>(2, {0,1,0,1 }) }));
	std::cout << dec2 << "\n";
	print_half_idempotent_relations<rational>(3);
	return 0;
}
