#include "HelperFunctions.hpp"
#include "TypeDefinitions.hpp"
#include "ModularArithmetic.hpp"

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

}