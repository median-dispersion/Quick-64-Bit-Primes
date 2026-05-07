#include "TypeDefinitions.hpp"
#include <iostream>
#include <chrono>
#include "Quick64BitPrimes.hpp"
#include <vector>

int main() {

    // Number to test
    q64bp::ui64 number;

    // Print a message
    std::cout << "Enter a number to test: ";

    // Get the number to test
    std::cin >> number;

    // Capture the test start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Check if the number is a prime using the Miller-Rabin primality test
    if (q64bp::millerRabinPrimalityTest(number)) {

        // Capture the test stop time
        auto stopTime = std::chrono::high_resolution_clock::now();

        // Calculate the duration between the start and stop time
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);

        // Print a message
        std::cout << number << " is prime!" << std::endl;

        // Print the prime factor
        std::cout << "Prime factors of " << number << ": " << number << "^1" << std::endl;

        // Print the total execution time
        std::cout << "Total execution time: " << duration.count() << " nanoseconds" << std::endl;

        // Exit
        return 0;

    }

    // Decompose the number into its prime factors using Pollard's rho algorithm
    std::vector<q64bp::PrimeFactor> primeFactors = q64bp::primeDecomposition(number);

    // Capture the test stop time
    auto stopTime = std::chrono::high_resolution_clock::now();

    // Calculate the duration between the start and stop time
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime);

    // Counter for counting the number of prime factors
    q64bp::ui64 counter = 0;

    // Print messages
    std::cout << number << " is not prime" << std::endl;
    std::cout << "Prime factors of " << number << ": ";

    // Loop through all prime factors
    for (auto& primeFactor : primeFactors) {

        // Print the prime factor
        std::cout << primeFactor.base << "^" << primeFactor.exponent;

        // Check if not the last prime factor
        if (counter < primeFactors.size() - 1) {

            // Print the separator
            std::cout << ", ";

        }

        // Increase the counter
        counter++;

    }

    // Print the total execution time
    std::cout << std::endl << "Total execution time: " << duration.count() << " nanoseconds" << std::endl;

    // Exit
    return 0;

}