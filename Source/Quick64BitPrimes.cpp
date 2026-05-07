#include "Quick64BitPrimes.hpp"
#include "TypeDefinitions.hpp"
#include "ModularArithmetic.hpp"
#include "HelperFunctions.hpp"
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>

namespace q64bp {

    // ============================================================================================
    // Check if a number is prime using the Miller-Rabin primality test
    // Based on: https://cp-algorithms.com/algebra/primality_tests.html#miller-rabin-primality-test
    // ============================================================================================
    bool millerRabinPrimalityTest(ui64 number) {

        // No number less than two can be prime
        if (number < 2) { return false; }

        // Check small primes directly witch is faster than using the Miller-Rabin primality test
        // Values from https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
        // These value are technically intended to be used as bases in the Miller-Rabin primality test
        // But they also happen to work great as a quick check for small primes
        for (ui64 prime : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {

            // If the number itself is the prime return true
            if (number == prime) { return true; }

            // If the number is divisible by the prime return false
            if (number % prime == 0) { return false; }

        }

        // Initialize the factor and exponent (d & s)
        ui64 factor = number - 1;
        ui64 exponent = 0;

        // Loop as long as factor is even using a bitwise AND check
        while (!(factor & 1)) {

            // Divide the factor by 2 using a right bit shift
            factor >>= 1;

            // Increase exponent
            exponent++;

        }

        // Deterministic set of bases that works for all 64-bit integers
        // Values from https://miller-rabin.appspot.com/
        for (ui64 base : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {

            // If the base is a multiple of the number continue with the next base
            if (base % number == 0) { continue; }

            // Calculate a result using modular exponentiation
            ui64 result = ModularArithmetic::exponentiation(base, factor, number);

            // If the result is 1 or number - 1 continue with the next base
            if (result == 1 || result == number - 1) { continue; }

            // Loop up to exponent - 1 times
            // Or until the result == number - 1
            for (ui64 loop = 1; loop < exponent && result != number - 1; loop++) {

                // Square the result using modular multiplication
                result = ModularArithmetic::multiplication(result, result, number);

            }

            // If the result is not number - 1 then the number is a composite so return false
            if (result != number - 1) { return false; }

        }

        // Return that the number is a prime
        return true;

    }

    // ============================================================================================
    // Prime decomposition using trial division
    // ============================================================================================
    void trialDivision(
        ui64 number,
        std::vector<ui64>& primes
    ) {

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // This is guaranteed by calling this function through primeDecomposition()

        // Loop as long as the number is divisible by two using a bitwise AND check
        while (!(number & 1)) {

            // Divide out the prime factor
            number /= 2;

            // Add the prime factor to the list of primes
            primes.push_back(2);

        }

        // Loop through every odd factor up to the square root of the number
        for(ui64 factor = 3; factor * factor <= number; factor += 2) {

            // Loop as long as the number is divisible by the factor
            while (number % factor == 0) {

                // Divide out the prime factor
                number /= factor;

                // Add the prime factor to the list of primes
                primes.push_back(factor);

            }

        }

        // If the number is more than one, it is a prime so add it to the list of primes
        if (number > 1) { primes.push_back(number); }

    }

    // ============================================================================================
    // Integer factorization using Pollard's rho algorithm and Floyd's cycle detection method
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Algorithm
    // ============================================================================================
    /* ui64 pollardFloydFactorization(ui64 number) {

        // This function is a legacy drop in replacement for pollardBrentFactorization()
        // It uses Floyd's cycle detection method and is therefore slightly slower

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // This is guaranteed by calling this function through primeDecomposition()

        // If the number is divisible by two return a factor of two
        if (!(number & 1)) { return 2; }

        // Initialize the random number generator
        static std::mt19937_64 rng(std::random_device{}());

        // Get a uniform distribution to choose a random starting postion for the variable of the polynomial function
        // The range is 2 to number - 2
        // This is to avoid trivial or degenerate cases
        std::uniform_int_distribution<ui64> variableDistribution(2, number - 2);

        // Get a uniform distribution to choose a random constant for the polynomial function
        // The range is 1 to number - 1
        // Defined as 1 to number - 2 because number - 2 will later be skipped by adding one
        std::uniform_int_distribution<ui64> constantDistribution(1, number - 2);

        // Loop until a nontrivial factor is found
        while (true) {

            // Set the tortoise to a random starting position in the allowed range
            ui64 tortoise = variableDistribution(rng);

            // Set the hare to the same starting position as the tortoise
            ui64 hare = tortoise;

            // Get a random value in the allowed range for the constant in the polynomial function
            ui64 constant = constantDistribution(rng);

            // Skip constant = number - 2
            // By adding one to the constant
            if (constant == number - 2) { constant++; }

            // Initialize the factor
            ui64 factor = 1;

            // Loop until a factor is found
            while (factor == 1) {

                // Advance the tortoise
                tortoise = HelperFunctions::modularPolynomial(tortoise, constant, number);

                // Advance the hare twice as much as the tortoise
                hare = HelperFunctions::modularPolynomial(HelperFunctions::modularPolynomial(hare, constant, number), constant, number);

                // Get the greatest common divisor between |tortoise - hare| and the number
                factor = std::gcd(tortoise < hare ? hare - tortoise : tortoise - hare, number);

            }

            // Return the factor if it is nontrivial
            if (factor < number) { return factor; }

            // If the factor is the number itself retry from the top with new random values for the polynomial function

        }

    } */

    // ============================================================================================
    // Integer factorization using Pollard's rho algorithm and Brent's cycle detection method
    // ============================================================================================
    ui64 pollardBrentFactorization(ui64 number) {

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // This is guaranteed by calling this function through primeDecomposition()

        // If the number is divisible by two return a factor of two
        if (!(number & 1)) { return 2; }

        // Initialize the random number generator
        static std::mt19937_64 rng(std::random_device{}());

        // Get a uniform distribution to choose a random starting postion for the variable of the polynomial function
        // The range is 2 to number - 2
        // This is to avoid trivial or degenerate cases
        std::uniform_int_distribution<ui64> variableDistribution(2, number - 2);

        // Get a uniform distribution to choose a random constant for the polynomial function
        // The range is 1 to number - 1
        // Defined as 1 to number - 2 because number - 2 will later be skipped by adding one
        std::uniform_int_distribution<ui64> constantDistribution(1, number - 2);

        // Loop until a nontrivial factor is found
        while (true) {

            // Set the tortoise to a random starting position in the allowed range
            ui64 tortoise = variableDistribution(rng);

            // Set the hare to the same starting position as the tortoise
            ui64 hare = tortoise;

            // Set the backup position to the position of the tortoise and hare
            ui64 backup = tortoise;

            // Get a random value in the allowed range for the constant in the polynomial function
            ui64 constant = constantDistribution(rng);

            // Skip constant = number - 2
            // By adding one to the constant
            if (constant == number - 2) { constant++; }

            // Batch size for Brent's optimization (controls how often the greatest common divisor is computed)
            // A value of 128 seems like a good value for all ranges
            ui64 batchSize = 128;

            // Initialize cycle length (doubles each phase)
            ui64 phaseLength = 1;

            // Initialize the product of differences (used for batched greatest common divisor calculations)
            ui64 product = 1;

            // Initialize the factor
            ui64 factor = 1;

            // Loop until a factor is found
            while (factor == 1) {

                // Teleport the tortoise to current hare position (start of new phase)
                tortoise = hare;

                // Advance hare by the phase length (explore sequence)
                for (ui64 index = 0; index < phaseLength; index++) {

                    // Advance the hare position
                    hare = HelperFunctions::modularPolynomial(hare, constant, number);

                }

                // Track how many steps have been taken
                ui64 step = 0;

                // Process the phase in batches
                while (step < phaseLength && factor == 1) {

                    // Save current hare position in case the fallback is required
                    backup = hare;

                    // Process up to the batch size or the remaining steps in the phase
                    for (ui64 index = 0; index < batchSize && index < phaseLength - step; index++) {

                        // Advance the hare position
                        hare = HelperFunctions::modularPolynomial(hare, constant, number);

                        // Multiply accumulated product by |tortoise - hare| mod number
                        product = ModularArithmetic::multiplication(
                            product,
                            tortoise < hare ? hare - tortoise : tortoise - hare,
                            number
                        );

                    }

                    // Compute the greatest common divisor of the accumulated product and the number
                    factor = std::gcd(product, number);

                    // Move forward by the batch size
                    step += batchSize;

                }

                // Double the phase length for next iteration using a left bit shift
                phaseLength <<= 1;

            }

            // Check if the found factor is the number itself
            // Do the fallback check
            if (factor == number) {

                do {

                    // Advance the backup position one step at a time
                    backup = HelperFunctions::modularPolynomial(backup, constant, number);

                    // Compute greatest common divisor of |tortoise - backup| and the number
                    factor = std::gcd(
                        tortoise < backup ? backup - tortoise : tortoise - backup,
                        number
                    );

                // Loop until a factor is found
                } while (factor == 1);

            }

            // Return the factor if it is nontrivial
            if (factor < number) { return factor; }

            // If the factor is still the number itself retry from the top with new random values for the polynomial function

        }

    }

    // ============================================================================================
    // Decompose a number into its prime factors
    // ============================================================================================
    std::vector<PrimeFactor> primeDecomposition(ui64 number) {

        // No number less than two can be broken down into prime factors
        if (number < 2) { return {}; }

        // Initialize a list for all factors and all primes
        std::vector<ui64> factors;
        std::vector<ui64> primes;

        // Reserve size for 64 entries
        // An unsigned 64-bit integer can only have 64 factors
        // 2⁶⁴ = 2 x 2 x 2 ... 64 times
        factors.reserve(64);
        primes.reserve(64);

        // Add the number to the list of factors
        factors.push_back(number);

        // Loop until there are no more factors
        while (!factors.empty()) {

            // Get the latest factor from the list of all factors
            ui64 factor = factors.back();

            // Remove the latest factor from the list of all factors
            factors.pop_back();

            // Check if the latest factor is a prime
            if (millerRabinPrimalityTest(factor)) {

                // Add the latest factor to the list of primes
                primes.push_back(factor);

                // Continue with the next factor
                continue;

            }

            // Check if the factor is less than 30'000'000
            if (factor < 30'000'000) {

                // Use trial division to break down the number into prime factors
                // Trial division is generally faster for small numbers up to around 30 million
                trialDivision(factor, primes);

                // Continue with the next factor
                continue;

            }

            // Break down the latest factor into a new factor and deduce the other new factor
            ui64 newFactor1 = pollardBrentFactorization(factor);
            ui64 newFactor2 = factor / newFactor1;

            // Add the two new factors to the list of all factors
            factors.push_back(newFactor1);
            factors.push_back(newFactor2);

        }

        // Sort the list of primes
        std::sort(primes.begin(), primes.end());

        // Initialize a list of all prime factors
        std::vector<PrimeFactor> primeFactors;

        // Reserve space for 15 prime factors
        // An unsigned 64-bit integer can only have 15 prime factors
        // Because the first 16 primes multiplied together would overflow the 64-bit range
        primeFactors.reserve(15);

        // Loop through the list of primes
        for (ui64 index = 0; index < primes.size();) {

            // Initialize the prime factor base and exponent
            ui64 base = primes[index];
            ui64 exponent = 0;

            // Check if the next prime is the same as the current prime
            while (base == primes[index] && index < primes.size()) {

                // Increase the exponent
                exponent++;

                // Continue with the next prime in the list
                index++;

            }

            // Add the prime factor
            primeFactors.emplace_back(base, exponent);

        }

        // Return the list of prime factors
        return primeFactors;

    }

    // ============================================================================================
    // Get the square root of a number modulo a prime using the Tonelli-Shanks algorithm
    // Based on the GO implementation: https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Go
    // ============================================================================================
    std::vector<ui64> tonelliShanksSquareRoot(
        ui64 number,
        ui64 prime
    ) {

        // If the provided value is not prime return no valid solutions
        if (!millerRabinPrimalityTest(prime)) { return {}; }

        // If the prime is two return the trivial solution
        if (prime == 2) { return {number % 2}; }

        // If the number exceeds the prime reduce it under modular arithmetic
        if (number >= prime) { number %= prime; }

        // If the number is zero return the trivial solution
        if (number == 0) { return {0}; }

        // If the number is a quadratic non residue return no solution
        if (ModularArithmetic::exponentiation(number, (prime - 1) >> 1, prime) != 1) { return {}; }

        // If the prime is in the form prime ≡ 3 (mod 4) use the fast path
        if ((prime & 3) == 3) {

            // Calculate the square root directly
            ui64 squareRoot = ModularArithmetic::exponentiation(number, (prime + 1) >> 2, prime);

            // Return the result
            return {squareRoot, prime - squareRoot};

        }

        // If the prime is in the form prime ≡ 5 (mod 8) use the fast path
        if ((prime & 7) == 5) {

            // For the case that prime ≡ 5 (mod 8) there are two different direct calculations
            // One for the case that the Legendre symbol is 1 and one for if it is -1
            // This path always assumes the Legendre symbol is 1, then calculates the first square root
            // It checks if the square root is valid and if not corrects the result to match the second case
            // After that it returns the valid square roots

            // Calculate the square root directly
            ui64 squareRoot = ModularArithmetic::exponentiation(number, (prime + 3) >> 3, prime);

            // Check if the square root is in the incorrect form
            if (ModularArithmetic::multiplication(squareRoot, squareRoot, prime) != number) {

                // Calculate the correction factor
                ui64 correction = ModularArithmetic::exponentiation(2, (prime - 1) >> 2, prime);

                // Apply the correction factor
                squareRoot = ModularArithmetic::multiplication(squareRoot, correction, prime);

            }

            // Return the result
            return {squareRoot, prime - squareRoot};

        }

        // Initialize the factor and exponent
        ui64 factor = prime - 1;
        ui64 exponent = 0;

        // Loop as long as factor is even using a bitwise AND check
        while (!(factor & 1)) {

            // Divide the factor by 2 using a right bit shift
            factor >>= 1;

            // Increase exponent
            exponent++;

        }

        // Initialize the quadratic non-residue
        ui64 quadraticNonResidue = 2;

        // Find a quadratic non-residue
        while (ModularArithmetic::exponentiation(quadraticNonResidue, (prime - 1) >> 1, prime) != prime - 1) { quadraticNonResidue++; }

        // Initialize the variables for the Tonelli-Shanks iteration
        ui64 squareRoot = ModularArithmetic::exponentiation(number, (factor + 1) >> 1, prime);
        ui64 currentFactor = ModularArithmetic::exponentiation(quadraticNonResidue, factor, prime);
        ui64 currentResidue = ModularArithmetic::exponentiation(number, factor, prime);
        ui64 currentExponent = exponent;

        // Loop until a solution is found
        while (currentResidue != 1) {

            // Initialize variables for exponent search
            ui64 temporaryFactor = currentResidue;
            ui64 newExponent = 0;

            // Find the smallest new exponent such that t^(2^i) ≡ 1 (mod p)
            while (temporaryFactor != 1 && newExponent + 1 < currentExponent) {

                temporaryFactor = ModularArithmetic::multiplication(temporaryFactor, temporaryFactor, prime);
                newExponent++;

            }

            // Initialize variables for factor search
            ui64 temporaryExponent = currentExponent - newExponent - 1;
            ui64 newFactor = currentFactor;

            // Find the new factor such that b = c^(2^(m - i - 1)) mod p
            while (temporaryExponent) {

                newFactor = ModularArithmetic::multiplication(newFactor, newFactor, prime);
                temporaryExponent--;

            }

            // Update variables according to the algorithm
            squareRoot = ModularArithmetic::multiplication(squareRoot, newFactor, prime);
            currentFactor = ModularArithmetic::multiplication(newFactor, newFactor, prime);
            currentResidue = ModularArithmetic::multiplication(currentResidue, currentFactor, prime);
            currentExponent = newExponent;

        }

        // Return the result
        return {squareRoot, prime - squareRoot};

    }

}