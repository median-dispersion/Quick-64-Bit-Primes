#include "Quick64BitPrimes/HelperFunctions.hpp"
#include "Quick64BitPrimes/TypeDefinitions.hpp"
#include "Quick64BitPrimes/ModularArithmetic.hpp"

namespace q64bp::HelperFunctions {

    // ============================================================================================
    // Modular polynomial f(x) = (x² + c) mod m
    // ============================================================================================
    ui64 modularPolynomial(
        ui64 variable,
        ui64 constant,
        ui64 modulus
    ) {
        return ModularArithmetic::addition(ModularArithmetic::multiplication(variable, variable, modulus), constant, modulus);
    }

    // ============================================================================================
    // Integer square root floor(sqrt(n))
    // Based on: https://en.wikipedia.org/wiki/Integer_square_root#Using_bitwise_operations
    // ============================================================================================
    ui64 integerSquareRoot(ui64 number) {

        // Initialize the square root
        ui64 squareRoot = 0;

        // Start with the highest power of four <= 2⁶⁴
        // 4611686018427387904 == static_cast<std::uint64_t>(1) << 62
        ui64 bit = 4611686018427387904;

        // Reduce the bit until it is the largest power of four <= number.
        while (bit > number) { bit >>= 2; }

        // Loop until all bits have been processed
        while (bit) {

            // Test whether adding this bit would keep the partial square root <= number.
            if (number >= squareRoot + bit) {

                // Subtract the trial value from the remainder
                number -= squareRoot + bit;

                // Shift current result down for next stage and add the bit
                squareRoot = (squareRoot >> 1) + bit;

            } else {

                // Reject the bit and just shift for next iteration
                squareRoot >>= 1;

            }

            // Move to the next lower power of four
            bit >>= 2;

        }

        // Return the square root
        return squareRoot;

    }

}