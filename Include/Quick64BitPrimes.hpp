#ifndef QUICK_64_BIT_PRIMES_HPP
#define QUICK_64_BIT_PRIMES_HPP

#include <cstdint>
#include <vector>

namespace q64bp {

    // ============================================================================================
    // Prime factor structure
    // ============================================================================================
    struct PrimeFactor {
        std::uint_fast64_t base;
        std::uint_fast64_t exponent;
    };

    // ============================================================================================
    // Check if a number is prime using the Miller-Rabin primality test
    // ============================================================================================
    bool millerRabinPrimalityTest(std::uint_fast64_t number);

    // ============================================================================================
    // Decompose a number into its prime factors
    // ============================================================================================
    std::vector<PrimeFactor> primeDecomposition(std::uint_fast64_t number);

    // ============================================================================================
    // Get the square root of a number modulo a prime using the Tonelli-Shanks algorithm
    // ============================================================================================
    std::vector<std::uint_fast64_t> tonelliShanksSquareRoot(
        std::uint_fast64_t number,
        std::uint_fast64_t prime
    );

}

#endif