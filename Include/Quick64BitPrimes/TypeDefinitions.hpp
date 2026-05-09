#ifndef QUICK_64_BIT_PRIMES_TYPE_DEFINITIONS_HPP
#define QUICK_64_BIT_PRIMES_TYPE_DEFINITIONS_HPP

#include <cstdint>

namespace q64bp {

    // Unsigned 64-bit integer alias
    using ui64 = std::uint_fast64_t;

    // Check if the __uint128_t type is available and set its alias
    #ifdef __SIZEOF_INT128__
    #define QUICK_64_BIT_PRIMES_UI128_AVAILABLE
    using ui128 = __uint128_t;
    #endif

}

#endif