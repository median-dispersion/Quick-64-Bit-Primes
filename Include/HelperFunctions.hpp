#ifndef QUICK_64_BIT_PRIMES_HELPER_FUNCTIONS_HPP
#define QUICK_64_BIT_PRIMES_HELPER_FUNCTIONS_HPP

#include "TypeDefinitions.hpp"

namespace q64bp::HelperFunctions {

    // ============================================================================================
    // Modular polynomial f(x) = (x² + c) mod m
    // ============================================================================================
    ui64 modularPolynomial(
        ui64 variable,
        ui64 constant,
        ui64 modulus
    );

}

#endif