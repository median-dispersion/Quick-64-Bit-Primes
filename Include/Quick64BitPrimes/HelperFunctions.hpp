#ifndef QUICK_64_BIT_PRIMES_HELPER_FUNCTIONS_HPP
#define QUICK_64_BIT_PRIMES_HELPER_FUNCTIONS_HPP

#include "Quick64BitPrimes/TypeDefinitions.hpp"

namespace q64bp::HelperFunctions {

    // ============================================================================================
    // Modular polynomial f(x) = (x² + c) mod m
    // ============================================================================================
    ui64 modularPolynomial(
        ui64 variable,
        ui64 constant,
        ui64 modulus
    );

    // ============================================================================================
    // Integer square root floor(sqrt(n))
    // ============================================================================================
    ui64 integerSquareRoot(ui64 number);

}

#endif