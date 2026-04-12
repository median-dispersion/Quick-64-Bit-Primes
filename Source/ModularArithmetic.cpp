#include "ModularArithmetic.hpp"
#include <cstdint>

namespace q64bp::ModularArithmetic {

    // ============================================================================================
    // Modular addition
    // ============================================================================================
    std::uint_fast64_t addition(
        std::uint_fast64_t addend1,
        std::uint_fast64_t addend2,
        std::uint_fast64_t modulus
    ) {
        return (static_cast<__uint128_t>(addend1) + addend2) % modulus;
    }

    // ============================================================================================
    // Modular multiplication
    // ============================================================================================
    std::uint_fast64_t multiplication(
        std::uint_fast64_t factor1,
        std::uint_fast64_t factor2,
        std::uint_fast64_t modulus
    ) {
        return static_cast<__uint128_t>(factor1) * factor2 % modulus;
    }

    // ============================================================================================
    // Modular exponentiation
    // https://en.wikipedia.org/wiki/Modular_exponentiation#Pseudocode
    // ============================================================================================
    std::uint_fast64_t exponentiation(
        std::uint_fast64_t base,
        std::uint_fast64_t exponent,
        std::uint_fast64_t modulus
    ) {

        // Initialize the result
        std::uint_fast64_t result = 1;

        // Reduce base to keep values small under modular arithmetic
        base %= modulus;

        // Loop until the exponent becomes zero
        while (exponent) {

            // Check if the exponent is odd using a bitwise AND
            if (exponent & 1) {

                // Multiply the result by the base mod modulus
                result = multiplication(result, base, modulus);

            }

            // Divide the exponent by 2 using a right bit shift
            exponent >>= 1;

            // Square the base mod modulus
            base = multiplication(base, base, modulus);

        }

        // Return the result
        return result;

    }

}