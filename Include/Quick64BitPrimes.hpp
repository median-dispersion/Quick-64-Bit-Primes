#ifndef QUICK_64_BIT_PRIMES_H
#define QUICK_64_BIT_PRIMES_H

#include <cstdint>

namespace q64bp {

    // Miller-Rabin primality test
    bool millerRabinPrimalityTest(std::uint_fast64_t number);

}

#endif