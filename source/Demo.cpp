#include "Half_Idempotent.h"
#include <random>
///@file
///@brief A demonstration file that can be compiled

using namespace Symmetric_Polynomials;


///Writes the twisted Pontryagin/symplectic classes pi_{s,j}/k_{s,j} in terms of the Chern classes under the forgetful map BU(n)->BSO(n) / hermitianization BU(n) -> BSp(n)
template<typename scalar_t>
void write_pontryagin_C2_in_terms_of_Chern_classes(int n)
{
	Half_Idempotent_Basis<scalar_t, Half_Idempotent_Variables<unsigned short>> hib(n);
	for (int s = 1; s <= n; s++)
		for (int i = 1; i <= n - s; i++) {
			auto twistedChern = hib.TwistedChern[s - 1][i - 1];
			Polynomial<scalar_t, Half_Idempotent_Variables<unsigned short>> twistedPontryagin(twistedChern.get_number_of_variables());
			for (auto it = twistedChern.begin(); it != twistedChern.end(); ++it) 
				twistedPontryagin.insert(it.exponent()+it.exponent(), it.coeff());
			std::cout << "k_{" << s << "," << i << "}= ";
			std::cout << hib(twistedPontryagin) << "\n";
		}
}


int main() {
	using namespace Symmetric_Polynomials;


	typedef unsigned char T; //change to short etc to allow for larger power


	std::cout << "The ring of symmetric polynomials on variables x_1,...,x_n,y_1,...,y_n with relations y_i^2=y_i can be minimally generated by the elements: \n";
	std::cout << "- The sum of idempotents a = y_1 + ... + y_n \n";
	std::cout << "- The Chern classes c_i that are the elementary symmetric polynomials on the x_1,...,x_n\n";
	std::cout << "- The twisted Chern classes c_{s,j}; each c_{s,j} is defined as the sum of all elements in the orbit of x_1....x_sy_{s+1}....y_{s+j} under the Sigma_n action\n";
	std::cout << "But the three types of classes satisfy relations. Enter the number n>=1 of variables x_1,...,x_n,y_1,...,y_n and all relations for the given n will be printed" << "\n";
	int n;
	std::cin >> n;
	if (n <= 0)
		std::cout << "Invalid n";
	else {
		std::cout << "The relations for n= " << n << " follow:" << "\n";
		print_half_idempotent_relations<Rational, Half_Idempotent_Variables<T>>(n,1);
	}

	//Write the Pontryagin/symplectic in terms of the Chern classes
	//write_pontryagin_C2_in_terms_of_Chern_classes<Rational>(10);

	//If we know the number of variables by compile time we can use that information with half_idempotent:
	//print_half_idempotent_relations<Rational, Half_Idempotent_Variables<T,4>>(2,1);
	//print_half_idempotent_relations<Rational, Half_Idempotent_Variables<T, 24>>(12);


	return 0;

}
