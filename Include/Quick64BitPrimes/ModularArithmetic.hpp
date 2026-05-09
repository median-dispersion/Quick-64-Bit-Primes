#ifndef QUICK_64_BIT_PRIMES_MODULAR_ARITHMETIC_HPP
#define QUICK_64_BIT_PRIMES_MODULAR_ARITHMETIC_HPP

#include "Quick64BitPrimes/TypeDefinitions.hpp"

namespace q64bp::ModularArithmetic {

    // ============================================================================================
    // Modular addition
    // ============================================================================================
    ui64 addition(
        ui64 addend1,
        ui64 addend2,
        ui64 modulus
    );

    // ============================================================================================
    // Modular multiplication
    // ============================================================================================
    ui64 multiplication(
        ui64 factor1,
        ui64 factor2,
        ui64 modulus
    );

    // ============================================================================================
    // Modular exponentiation
    // ============================================================================================
    ui64 exponentiation(
        ui64 base,
        ui64 exponent,
        ui64 modulus
    );

}

#endif