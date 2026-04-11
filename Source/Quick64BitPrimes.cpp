#include "Quick64BitPrimes.hpp"
#include <initializer_list>
#include "ModularArithmetic.hpp"

namespace q64bp {

    // ============================================================================================
    // Miller-Rabin primality test
    // https://cp-algorithms.com/algebra/primality_tests.html#miller-rabin-primality-test
    // ============================================================================================
    bool millerRabinPrimalityTest(std::uint_fast64_t number) {

        // Check if the number is less than 2
        if (number < 2) {

            // No number less than 2 can be prime
            // Return false
            return false;

        }

        // Check small primes directly witch is faster than using the Miller-Rabin primality test
        // Values from https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
        // These value are technically intended to be used as bases in the Miller-Rabin primality test
        // But they also happen to work great as a quick check for small primes
        for (std::uint_fast64_t prime : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {

            // Check if the number is divisible by the prime
            if (number % prime == 0) {

                // Return true only if the number itself is the prime else return false
                return number == prime;

            }

        }

        // Initialize the factor and exponent
        std::uint_fast64_t factor = number - 1;
        std::uint_fast64_t exponent = 0;

        // Loop as long as factor is even using a bitwise AND check
        while ((factor & 1) == 0) {

            // Divide the factor by 2 using a right bit shift
            factor >>= 1;

            // Increase exponent
            exponent++;

        }

        // Deterministic set of bases that works for all 64-bit integers
        // Values from https://miller-rabin.appspot.com/
        for (std::uint_fast64_t base : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {

            // Check if the base is a multiple of the number
            if (base % number == 0) {

                // Skip this base
                continue;

            }

            // Calculate a result using modular exponentiation
            std::uint_fast64_t result = ModularArithmetic::exponentiation(base, factor, number);

            // Check if the result is 1 or number - 1
            if (result == 1 || result == number - 1) {

                // Skip this base
                continue;

            }

            // Flag for checking if the number is a composite
            bool composite = true;

            // Repeat squaring up to exponent - 1 times
            for (std::uint_fast64_t repeat = 1; repeat < exponent; repeat++) {

                // Square the result using modular multiplication
                result = ModularArithmetic::multiplication(result, result, number);

                // Check if the result is number - 1
                if (result == number - 1) {

                    // Set the composite flag to false
                    composite = false;

                    // Stop the squaring loop immediately
                    break;

                }

            }

            // Check if the number is a composite
            if (composite) {

                // Return that the number is not a prime
                return false;

            }

        }

        // Return that the number is a prime
        return true;

    }

}