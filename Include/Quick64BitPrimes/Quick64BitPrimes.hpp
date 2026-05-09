#ifndef QUICK_64_BIT_PRIMES_HPP
#define QUICK_64_BIT_PRIMES_HPP

#include "Quick64BitPrimes/TypeDefinitions.hpp"
#include <vector>

namespace q64bp {

    // ============================================================================================
    // Prime factor structure
    // ============================================================================================
    struct PrimeFactor {
        ui64 base;
        ui64 exponent;
    };

    // ============================================================================================
    // Check if a number is prime using the Miller-Rabin primality test
    // ============================================================================================
    bool millerRabinPrimalityTest(ui64 number);

    // ============================================================================================
    // Decompose a number into its prime factors
    // ============================================================================================
    std::vector<PrimeFactor> primeDecomposition(ui64 number);

    // ============================================================================================
    // Get the square root of a number modulo a prime using the Tonelli-Shanks algorithm
    // ============================================================================================
    std::vector<ui64> tonelliShanksSquareRoot(
        ui64 number,
        ui64 prime
    );

    // ============================================================================================
    // Get Fermat's sum of two squares representation of a prime
    // ============================================================================================
    std::vector<ui64> fermatSumOfTwoSquares(ui64 prime);

}

#endif