#ifndef QUICK_64_BIT_PRIMES_HPP
#define QUICK_64_BIT_PRIMES_HPP

#include "Quick64BitPrimes/TypeDefinitions.hpp"
#include <vector>
#include <optional>
#include <utility>

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
    // Get the square roots of a number modulo a prime using the Tonelli-Shanks algorithm (r² ≡ n (mod p))
    // ============================================================================================
    std::optional<std::pair<ui64, std::optional<ui64>>> tonelliShanksAlgorithm(
        ui64 number,
        ui64 prime
    );

    // ============================================================================================
    // Get Fermat's sum of two squares representation of a prime (x² + y² = p)
    // ============================================================================================
    std::optional<std::pair<ui64, ui64>> fermatSumOfTwoSquaresTheorem(ui64 prime);

}

#endif