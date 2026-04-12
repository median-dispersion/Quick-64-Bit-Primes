#include "Quick64BitPrimes.hpp"
#include <cstdint>
#include <initializer_list>
#include "ModularArithmetic.hpp"
#include <random>
#include <numeric>
#include <vector>

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

    // ============================================================================================
    // Polynomial function f(x) = (x² + c) mod m
    // ============================================================================================
    std::uint_fast64_t polynomial(
        std::uint_fast64_t variable,
        std::uint_fast64_t constant,
        std::uint_fast64_t modulus
    ) {
        return ModularArithmetic::addition(ModularArithmetic::multiplication(variable, variable, modulus), constant, modulus);
    }

    // ============================================================================================
    // Integer factorization using Pollard's rho algorithm and Floyd's cycle detection method
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Algorithm
    // ============================================================================================
    std::uint_fast64_t pollardFloydFactorization(std::uint_fast64_t number) {

        // Do a quick safety check if the number is divisible by two
        if (number % 2 == 0) {

            // Return a factor of two
            return 2;

        }

        // Initialize the random number generator
        static std::mt19937_64 rng(std::random_device{}());

        // Get a uniform distribution to choose a random starting postion for the tortoise
        // The range is 2 to number - 2
        // This is to avoid trivial or degenerate cases
        std::uniform_int_distribution<std::uint_fast64_t> tortoiseDistribution(2, number - 2);

        // Get a uniform distribution to choose a random constant for the polynomial function
        // The range is 1 to number - 1
        std::uniform_int_distribution<std::uint_fast64_t> constantDistribution(1, number - 1);

        // Loop until a nontrivial factor is found
        while (true) {

            // Set the tortoise to a random starting position in the allowed range
            std::uint_fast64_t tortoise = tortoiseDistribution(rng);

            // Set the hare to the same starting position as the tortoise
            std::uint_fast64_t hare = tortoise;

            // Get a random value in the allowed range for the constant in the polynomial function
            std::uint_fast64_t constant = constantDistribution(rng);

            // Check if the constant is number - 2
            // Loop until it isn't
            // A constant of number - 2 can cause issues for example shot cycles and symmetry
            while (constant == number - 2) {

                // Choose a new random constant
                constant = constantDistribution(rng);

            }

            // Initialize the factor
            std::uint_fast64_t factor = 1;

            // Loop until a factor is found
            while (factor == 1) {

                // Advance the tortoise
                tortoise = polynomial(tortoise, constant, number);

                // Advance the hare twice as much as the tortoise
                hare = polynomial(polynomial(hare, constant, number), constant, number);

                // Get the greatest common divisor between |tortoise - hare| and the number
                factor = std::gcd((tortoise < hare) ? hare - tortoise : tortoise - hare, number);

            }

            // Check if the factor is not the number itself
            if (factor < number) {

                // Return the nontrivial factor
                return factor;

            }

            // If the factor is the number itself retry from the top with new random values for the polynomial function

        }

    }

    // ============================================================================================
    // Prime decomposition
    // ============================================================================================
    std::vector<PrimeFactor> primeDecomposition(std::uint_fast64_t number) {

        // Check if the number is less than 2
        if (number < 2) {

            // No number less that two can be a prime or broken down into prime factors
            // Return no prime factors
            return {};

        }

        // Check if the number itself is prime
        if (millerRabinPrimalityTest(number)) {

            // No prime can be broken down into prime factors
            // Return no prime factors
            return {};

        }

        // Initialize a list of all factors of the number
        std::vector<std::uint_fast64_t> allFactors;

        // Initialize a list of all prime factors of the number
        std::vector<PrimeFactor> primeFactors;

        // Get a factor of the number
        std::uint_fast64_t factor1 = pollardFloydFactorization(number);

        // Deduce the other factor
        std::uint_fast64_t factor2 = number / factor1;

        // Add both factors to the list of all factors
        allFactors.push_back(factor1);
        allFactors.push_back(factor2);

        // Loop until there are no more none prime factors
        while (!allFactors.empty()) {

            // Get the latest factor from the list of all factors
            std::uint_fast64_t latestFactor = allFactors.back();

            // Remove the latest factor from the list of all factors
            allFactors.pop_back();

            // Check if the latest factor is one
            if (latestFactor == 1) {

                // One is not a prime and can not be broken down
                // Continue with next factor
                continue;

            }

            // Check if the latest factor is a prime
            if (millerRabinPrimalityTest(latestFactor)) {

                // Initialize the prime factor exponent
                std::uint_fast64_t exponent = 1;

                // Loop through all other factors of the number
                for (auto& otherFactor : allFactors) {

                    // Loop until the other factor is no longer divisible by the prime factor
                    while (otherFactor % latestFactor == 0) {

                        // Divide out the prime factor
                        otherFactor /= latestFactor;

                        // Increase the prime factor exponent
                        exponent++;

                    }

                }

                // Add the prime factor to the list of all prime factors
                primeFactors.emplace_back(latestFactor, exponent);

                // Continue with the next factor
                continue;

            }

            // Break down the latest factor into a new factor
            std::uint_fast64_t newFactor1 = pollardFloydFactorization(latestFactor);

            // Deduce the other factor
            std::uint_fast64_t newFactor2 = latestFactor / newFactor1;

            // Add the two new factors into the list of all factors
            allFactors.push_back(newFactor1);
            allFactors.push_back(newFactor2);

        }

        // Return all prime factors
        return primeFactors;

    }

}