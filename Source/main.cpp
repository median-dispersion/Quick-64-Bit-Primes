#include <cstdint>
#include <iostream>
#include <chrono>
#include "Quick64BitPrimes.hpp"

int main() {

    // Number to test
    std::uint_fast64_t number;

    // Print a message
    std::cout << "Enter a number to test: ";

    // Get the number to test
    std::cin >> number;

    // Capture the test start time
    auto start = std::chrono::high_resolution_clock::now();

    // Check if the number is a prime using the Miller-Rabin primality test
    bool isPrime = q64bp::millerRabinPrimalityTest(number);

    // Capture the test stop time
    auto stop = std::chrono::high_resolution_clock::now();

    // Calculate the duration between the start and stop time
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    // If the number is prime
    if (isPrime) {

        // Print a message
        std::cout << number << " is a prime!" << std::endl;

    // If the number is not a prime
    } else {

        // Print a message
        std::cout << number << " is not a prime" << std::endl;

    }

    // Print the total test time
    std::cout << "Total test duration: " << duration.count() << " nanoseconds" << std::endl;

    // Exit
    return 0;

}