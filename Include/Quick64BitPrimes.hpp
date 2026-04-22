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
    // Miller-Rabin primality test
    // ============================================================================================
    bool millerRabinPrimalityTest(std::uint_fast64_t number);

    // ============================================================================================
    // Prime decomposition
    // ============================================================================================
    std::vector<PrimeFactor> primeDecomposition(std::uint_fast64_t number);

}

#endif