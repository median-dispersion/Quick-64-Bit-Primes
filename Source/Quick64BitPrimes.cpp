#include "Quick64BitPrimes/Quick64BitPrimes.hpp"
#include "Quick64BitPrimes/TypeDefinitions.hpp"
#include "Quick64BitPrimes/ModularArithmetic.hpp"
#include "Quick64BitPrimes/HelperFunctions.hpp"
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <optional>
#include <utility>

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
        // Write number - 1 as factor^exponent
        ui64 factor = number - 1;
        ui64 exponent = 0;

        // Loop as long as factor is even using a bitwise AND check
        // Factoring out all powers of two from number - 1
        while (!(factor & 1)) {

            // Divide the factor by 2 using a right bit shift
            factor >>= 1;

            // Increase exponent
            exponent++;

        }

        // Deterministic set of bases that works for all 64-bit integers
        // Values from https://miller-rabin.appspot.com/
        // If none of these bases prove the number is a composite the number must be prime
        for (ui64 base : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {

            // If the base is a multiple of the number continue with the next base
            // This is to avoid cases where the modular exponentiation would result in 0 instead of 1 or number - 1
            // In that case it would provide no useful information for detecting primality but also not proving the number composite
            if (base % number == 0) { continue; }

            // Calculate the starting point of the Miller-Rabin primality test
            ui64 result = ModularArithmetic::exponentiation(base, factor, number);

            // If the result is 1 or number - 1, continue with the next base
            // In the case that result is 1, squaring it repeatedly will make it stay 1
            // In the case that result is number - 1, squaring it would make it 1, witch will always stay 1
            // In those scenarios the base can't be used to prove the number is a composite
            // Indicating that it might be prime, so continue the test with the next base
            if (result == 1 || result == number - 1) { continue; }

            // Loop up to exponent - 1 times
            // Or until the result == number - 1
            for (ui64 loop = 1; loop < exponent && result != number - 1; loop++) {

                // Square the result using modular multiplication
                // If the result becomes number - 1, the base cant be used to prove the number is a composite
                // Indicating that it might be prime, so continue the test with the next base
                result = ModularArithmetic::multiplication(result, result, number);

            }

            // If the result is not number - 1, then the number is proven to be a composite, so return false
            // This is because it did not behave like a prime under modular arithmetic eventually becoming number - 1
            if (result != number - 1) { return false; }

        }

        // If no base could prove the number is a composite, the number must be prime
        // Return that the number is a prime
        return true;

    }

    // ============================================================================================
    // Prime decomposition using trial division
    // ============================================================================================
    void trialDivisionUnsafe(
        ui64 number,
        std::vector<ui64>& primes
    ) {

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // Otherwise it is unsafe to use!
        // These checks are guaranteed by calling this function through primeDecomposition()

        // Loop as long as the number is divisible by two using a bitwise AND check
        while (!(number & 1)) {

            // Divide out the prime factor using a right bit shift
            number >>= 1;

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
    /* ui64 pollardFloydFactorizationUnsafe(ui64 number) {

        // This function is a legacy drop in replacement for pollardBrentFactorizationUnsafe()
        // It uses Floyd's cycle detection method and is therefore slightly slower
        // It is also not hardened against infinite loops and could hang forever!
        // It shouldn't be used and is just kept as a reference

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // Otherwise it is unsafe to use!
        // These checks are guaranteed by calling this function through primeDecomposition()

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
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Algorithm
    // ============================================================================================
    ui64 pollardBrentFactorizationUnsafe(ui64 number) {

        // Any number passed to this function must be:
        // number >= 2
        // number != prime
        // Otherwise it is unsafe to use!
        // These checks are guaranteed by calling this function through primeDecomposition()

        // If that is the case the function is practically "guaranteed" to return a non trivial factor
        // However there still is an astronomically small change the algorithm will fail
        // In that case the function will return 0 witch is an invalid factor!
        // If that occurs a deterministic approach like trial division must be used to find the factors

        // If the number is divisible by two return a factor of two
        if (!(number & 1)) { return 2; }

        // Initialize the random number generator
        static std::mt19937_64 rng(std::random_device{}());

        // Get a uniform distribution to choose a random constant for the polynomial function
        // The range is 1 to number - 1
        // Defined as 1 to number - 2, because number - 2 will later be skipped by adding one
        std::uniform_int_distribution<ui64> anchor(1, number - 2);

        // Get a uniform distribution to choose a random start postion for the tortoise, hare and fox
        // The range is 2 to number - 2
        // This is to avoid trivial or degenerate cases
        std::uniform_int_distribution<ui64> position(2, number - 2);

        // Define the maximum number of steps the hare can do in one sprint
        // This is Brent's tuning factor "m"
        // A value of 128 seems like a good value for all input ranges
        ui64 sprintLength = 128;

        // Limit the number of attempts to find a non trivial factor
        // Pollard's Rho "should" find a non trivial factor
        // Especially when retrying with randomized parameters for the polynomial function
        // But there seems to be no proof that if all possible combination of parameters are exhausted that a non trivial factor will be found
        // This is why the number of attempts is limited and if no factor is found an invalid solution of 0 is returned
        // During testing the maximum number of attempts reached was 32
        // This multiplied with a safety factor of 8
        // So for 99.999..% of the time the algorithm will return a valid solution
        // And for the astronomically small chance that it fails it limits the number of attempts instead of hanging forever
        for (ui64 attempt = 0; attempt < 256; attempt++) {

            // Initialize the factor as one
            ui64 factor = 1;

            // Initialize the start positions of the hare, tortoise and fox
            // Randomize the start position for each attempt
            ui64 hare = position(rng);
            ui64 tortoise = hare;
            ui64 fox = hare;

            // Initialize the constant of the polynomial function
            // Randomize the constant for each attempt
            ui64 constant = anchor(rng);

            // Skip constant = number - 2, by adding one to the constant
            if (constant == number - 2) { constant++; }

            // Repeatedly perform races between the tortoise and the hare
            // Doubling the race length each time
            // Do this until the maximum race length has been reached or a factor has been found
            // During testing the maximum length the race reached was 262'144
            // This multiplied with a safety factor of 8
            for (ui64 raceLength = 1; raceLength < 2'097'152 && factor == 1; raceLength <<= 1) {

                // The tortoise and hare start on the same position
                tortoise = hare;

                // Let the hare take as many steps as the race is long (cheater!)
                for (ui64 step = 0; step < raceLength; step++) {

                    // Advance the position of the hare
                    hare = HelperFunctions::modularPolynomial(hare, constant, number);

                }

                // Initialize the accumulative distance between the tortoise and the hare
                // This will be updated after each sprint of the hare
                ui64 distance = 1;

                // Let the hare perform a short sprint
                // Do as many sprints as fit in the race or until a factor is found
                for (ui64 sprint = 0; sprint < raceLength && factor == 1; sprint += sprintLength) {

                    // Before each sprint, the fox and hare meat up
                    // The fox stays at this position, keeping a lookout for the tortoise
                    fox = hare;

                    // Let the hare take as steps many as the sprint is long
                    // Or as many steps as are remaining in the race
                    for (ui64 step = 0; step < sprintLength && step < raceLength - sprint; step++) {

                        // Advance the position of the hare
                        hare = HelperFunctions::modularPolynomial(hare, constant, number);

                        // Update the accumulative distance between the tortoise and hare
                        distance = ModularArithmetic::multiplication(
                            distance,
                            tortoise < hare ? hare - tortoise : tortoise - hare,
                            number
                        );

                    }

                    // After each sprint try to find a hidden factor in the accumulate distance
                    factor = std::gcd(distance, number);

                }

            }

            // After all the races are completed, check if the factor is the number itself
            // This means the hare overshot the goal of the race
            if (factor == number) {

                // Reset the factor to try and find a non trivial one
                factor = 1;

                // The fox is a good sport and performs one final race between him and the tortoise
                // He takes as many steps as this final race is long, or until a factor has been found
                // During testing the maximum length of the final race reached was 128
                // This multiplied with a safety factor of 8
                for (ui64 step = 0; step < 1024 && factor == 1; step++) {

                    // Advance the position of the fox
                    // The fox is only doing on step at a time instead of sprinting like tha hare did
                    fox = HelperFunctions::modularPolynomial(fox, constant, number);

                    // After each step try to find a hidden factor in the distance between the tortoise and fox
                    factor = std::gcd(
                        tortoise < fox ? fox - tortoise : tortoise - fox,
                        number
                    );

                }

            }

            // If a non trivial factor was found return it
            if (factor > 1 && factor < number) { return factor; }

            // If no non trivial factor was found
            // Retry with new randomized start positions and parameters for the polynomial function

        }

        // If the maximum number of attempts was reached (extremely unlikely)
        // No non trivial factor was found
        // Return an invalid solution
        return 0;

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
                trialDivisionUnsafe(factor, primes);

                // Continue with the next factor
                continue;

            }

            // Break down the latest factor into a new factor using Pollard's rho
            ui64 newFactor1 = pollardBrentFactorizationUnsafe(factor);

            // Check if a valid new factor was found
            if (newFactor1) {

                // Deduce the other new factor
                ui64 newFactor2 = factor / newFactor1;

                // Add the two new factors to the list of all factors
                factors.push_back(newFactor1);
                factors.push_back(newFactor2);

            // If Pollard's rho failed to find a non trivial factor
            } else {

                // Use a deterministic approach like trial division
                // This is significantly slower especially for large factors
                // But in practice this should never happen and is just a fallback
                trialDivisionUnsafe(factor, primes);

            }

        }

        // Sort the list of primes from smallest to largest
        std::sort(primes.begin(), primes.end());

        // Initialize a list of all prime factors
        std::vector<PrimeFactor> primeFactors;

        // Reserve space for 15 prime factors
        // An unsigned 64-bit integer can only have 15 prime factors
        // Because the first 16 primes multiplied together would overflow the 64-bit range
        primeFactors.reserve(15);

        // Loop through the list of sorted primes
        for (ui64 index = 0; index < primes.size();) {

            // Initialize the prime factor base and exponent
            ui64 base = primes[index];
            ui64 exponent = 0;

            // Check if the next prime is the same as the current prime
            while (base == primes[index] && index < primes.size()) {

                // Increase the exponent of the current prime
                exponent++;

                // Continue with the next prime in the list
                index++;

            }

            // Add the prime factor to the list of prime factors
            primeFactors.emplace_back(base, exponent);

        }

        // Return the list of prime factors
        return primeFactors;

    }

    // ============================================================================================
    // Get the square roots of a number modulo a prime using the Tonelli-Shanks algorithm (r² ≡ n (mod p))
    // Based on the GO implementation: https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Go
    // ============================================================================================
    std::pair<ui64, ui64> tonelliShanksAlgorithmUnsafe(
        ui64 number,
        ui64 prime
    ) {

        // Any input passed to this function must be:
        // number % prime != 0
        // prime == must be prime
        // prime == must be odd
        // Legendre symbol == 1
        // Otherwise it is unsafe to use!
        // These checks are guaranteed by calling this function through tonelliShanksAlgorithm()

        // If the prime is in the form prime ≡ 3 (mod 4) use the fast path
        if ((prime & 3) == 3) {

            // Calculate the square root directly
            // Prime + 1 is safe because the largest 64-bit prime + 1 will not overflow
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
            // Prime + 3 is safe because the largest 64-bit prime + 3 will not overflow
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
        // Writing prime - 1 as factor^exponent
        ui64 factor = prime - 1;
        ui64 exponent = 0;

        // Loop as long as factor is even using a bitwise AND check
        // Factoring out all powers of two from prime - 1
        while (!(factor & 1)) {

            // Divide the factor by 2 using a right bit shift
            factor >>= 1;

            // Increase exponent
            exponent++;

        }

        // Initialize the quadratic non-residue
        ui64 quadraticNonResidue = 2;

        // Find a quadratic non-residue by repeatedly checking if the Legendre symbol is not -1
        // Increasing the quadratic non-residue variable by one each loop and checking again
        // For an odd prime a quadratic non-residue must exist
        while (ModularArithmetic::exponentiation(quadraticNonResidue, (prime - 1) >> 1, prime) != prime - 1) { quadraticNonResidue++; }

        // Initialize the variables for the Tonelli-Shanks iteration
        ui64 squareRoot = ModularArithmetic::exponentiation(number, (factor + 1) >> 1, prime);
        ui64 currentFactor = ModularArithmetic::exponentiation(quadraticNonResidue, factor, prime);
        ui64 currentResidue = ModularArithmetic::exponentiation(number, factor, prime);
        ui64 currentExponent = exponent;

        // Loop until a solution is found
        // With a Legendre symbol of 1 and an odd prime a solution must exist so this will terminate
        while (currentResidue != 1) {

            // Initialize variables for exponent search
            ui64 temporaryFactor = currentResidue;
            ui64 newExponent = 0;

            // Find the smallest new exponent such that t^(2^i) ≡ 1 (mod p)
            while (temporaryFactor != 1 && newExponent + 1 < currentExponent) {

                temporaryFactor = ModularArithmetic::multiplication(temporaryFactor, temporaryFactor, prime);
                newExponent++;

            }

            // This should not be necessary as long as the inputs are validated
            // Meaning the Legendre symbol of the number == 1 and the prime being odd
            // if (!currentExponent) { throw std::runtime_error("The Tonelli-Shanks algorithm failed!"); }

            // Initialize variables for factor search
            // For valid Tonelli-Shanks inputs, currentExponent can never be 0, so no underflow or infinite loop risk
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

    // ============================================================================================
    // Get the square roots of a number modulo a prime using the Tonelli-Shanks algorithm (r² ≡ n (mod p))
    // Input validation before the main algorithm
    // ============================================================================================
    std::optional<std::pair<ui64, std::optional<ui64>>> tonelliShanksAlgorithm(
        ui64 number,
        ui64 prime
    ) {

        // If the provided value is not prime return no valid solutions
        if (!millerRabinPrimalityTest(prime)) { return {}; }

        // If the prime is two return the trivial solution
        if (prime == 2) { return {{number & 1, {}}}; }

        // If the number exceeds the prime reduce it under modular arithmetic
        if (number >= prime) { number %= prime; }

        // If the number is zero return the trivial solution
        if (number == 0) { return {{0, {}}}; }

        // If the Legendre symbol is not one return no solution
        if (ModularArithmetic::exponentiation(number, (prime - 1) >> 1, prime) != 1) { return {}; }

        // Run the main Tonelli-Shanks algorithm and return the result
        return tonelliShanksAlgorithmUnsafe(number, prime);

    }

    // ============================================================================================
    // Get Fermat's sum of two squares representation of a prime (x² + y² = p)
    // ============================================================================================
    std::optional<std::pair<ui64, ui64>> fermatSumOfTwoSquaresTheorem(ui64 prime) {

        // If the provided value is not prime return no valid solutions
        if (!millerRabinPrimalityTest(prime)) { return {}; }

        // If the prime is two return the trivial solutions
        if (prime == 2) { return {{1, 1}}; }

        // If the prime in not in the form prime ≡ 1 (mod 4) return no valid solutions
        if ((prime & 3) != 1) { return {}; }

        // Set x to the prime
        ui64 x = prime;

        // Get r² ≡ p-1 (mod p) using Tonelli-Shanks
        // Use the "unsafe" function because at this point the prime is validated to be safe
        // Set y to the first result of the Tonelli-Shanks algorithm
        // For an odd prime in the form prime ≡ 1 (mod 4) at least one result is guaranteed
        ui64 y = tonelliShanksAlgorithmUnsafe(prime - 1, prime).first;

        // Run the Euclidean algorithm until y <= √prime
        while (y > prime / y) {

            // Get the remainder of x divided by y
            ui64 remainder = x % y;

            // Set x to the current divisor
            x = y;

            // Set y to the remainder
            y = remainder;

        }

        // Return the result
        return {{HelperFunctions::integerSquareRoot(prime - y * y), y}};

    }

}