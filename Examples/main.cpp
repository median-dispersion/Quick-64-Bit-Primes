#include "Quick64BitPrimes/Quick64BitPrimes.hpp"
#include <iostream>
#include <chrono>
#include <vector>

int main() {

    // Get the input number
    q64bp::ui64 number;
    std::cout << "Enter a number to test: ";
    std::cin >> number;

    // Capture the test start time
    auto start_time = std::chrono::high_resolution_clock::now();

    // Check if the number is a prime using the Miller-Rabin primality test
    if (q64bp::miller_rabin_primality_test(number)) {

        // Capture the total execution time
        auto stop_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time);

        // Print messages
        std::cout << number << " is prime!" << std::endl;
        std::cout << "Prime factors of " << number << ": " << number << "^1" << std::endl;
        std::cout << "Total execution time: " << duration.count() << " nanoseconds" << std::endl;

        // Exit
        return 0;

    }

    // Decompose the number into its prime factors using Pollard's rho algorithm
    std::vector<q64bp::PrimeFactor> prime_factors = q64bp::prime_decomposition(number);

    // Capture the total execution time
    auto stop_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time);

    // Counter for counting the number of prime factors
    q64bp::ui64 counter = 0;

    // Print messages
    std::cout << number << " is not prime" << std::endl;
    std::cout << "Prime factors of " << number << ": ";

    // Loop through all prime factors
    for (auto& prime_factor : prime_factors) {

        // Print the prime factor
        std::cout << prime_factor.base << "^" << prime_factor.exponent;

        // Print the separator if the isn't the last prime factor
        counter++; if (counter < prime_factors.size()) { std::cout << ", "; }

    }

    // Print the total execution time
    std::cout << std::endl << "Total execution time: " << duration.count() << " nanoseconds" << std::endl;

    // Exit
    return 0;

}