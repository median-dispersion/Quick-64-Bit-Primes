# Quick 64-Bit Primes

<p align="center">
    <a href="./Assets/Banner.svg">
        <img src="./Assets/Banner.svg" style="width: 900px;">
    </a>
</p>

<a href="./LICENSE"><img alt="GitHub License" src="https://img.shields.io/github/license/median-dispersion/Quick-64-Bit-Primes?style=for-the-badge"></a>
<a href="https://github.com/median-dispersion/Quick-64-Bit-Primes/releases/latest"><img alt="GitHub Release" src="https://img.shields.io/github/v/release/median-dispersion/Quick-64-Bit-Primes?style=for-the-badge"></a>

This repository contains C++ implementations of the following prime number algorithms:

- [Miller-Rabin primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)
- [Pollard's rho algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm) using [Brent's cycle detection](https://en.wikipedia.org/wiki/Cycle_detection#Brent.27s_algorithm)
- [Tonelli-Shanks algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm)
- [Fermat's sum of two squares theorem](https://en.wikipedia.org/wiki/Fermat%27s_theorem_on_sums_of_two_squares)

All algorithm implementations are safe, deterministic and accept the full range of unsigned 64-bit integers as a number between $0$ and $2^{64}-1$. Depending on the algorithm, the runtime can range from just a few nanoseconds up to about 20 microseconds in the worst case for a full 64-bit value.

# Algorithms

## Miller-Rabin primality test

To check if a number is prime, the Miller-Rabin primality test can be used.

The Miller-Rabin primality test is a probabilistic algorithm that checks whether a number is prime by testing if it behaves like a prime under repeated modular exponentiation. If it fails any test, it’s guaranteed to be composite, but if it passes multiple rounds, it’s very likely prime. Its time complexity is $O(k\log^{3}n)$, where $n$ is the number being tested and $k$ is the number of rounds.

Probabilistic, in this context, means the test does not always produce a correct result and may occasionally incorrectly classify a number as prime. However, this particular implementation is both deterministic and fully accurate. Restricting the input to the unsigned 64-bit integer range allows the use of a fixed set of bases that guarantees 100% correct primality testing within that range. It also simplifies the complexity to $O(\log^{3}n)$.

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any number within the 64-bit integer range can be checked in less than 110 nanoseconds. The slowest cases occur near the upper end of the 64-bit range, where the test takes around 110 nanoseconds. The benchmark results can be found in [Documentation/Primality test data.csv](./Documentation/Primality%20test%20data.csv).

<p align="center">
    <a href="./Documentation/Primality test chart.svg">
        <img src="./Documentation/Primality test chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
bool q64bp::miller_rabin_primality_test(q64bp::ui64 number);
```

#### Example

```c++
#include "Quick64BitPrimes/miller_rabin_primality_test.hpp"
#include <iostream>

int main() {

    int number = 123;

    if (q64bp::miller_rabin_primality_test(number)) {

        std::cout << number << " is prime!" << std::endl;

    } else {

        std::cout << number << " is not prime." << std::endl;

    }

}
```

## Pollard's rho algorithm

To decompose a number into its prime factors, Brent's improved variant of Pollard's rho can be used.

Pollard's rho is a probabilistic algorithm that finds a nontrivial factor of a composite number $n$ by generating a pseudo-random sequence modulo $n$ and using cycle detection to reveal a hidden common divisor. Brent's variant of Pollard's rho is a cycle-detection improvement over Floyd's method that batches iterations to reduce the expensive greatest common divisor computations while still finding a nontrivial factor of a composite number. Its worst-case time complexity is $O(n^{1/4})$, and in practice it's often faster.

Probabilistic, in this context, means that the algorithm might fail to find a factor that is not one or the number itself. However, this particular implementation of Pollard’s rho behaves deterministically in practice. Whenever it fails, it simply retries with a new set of random values for the polynomial function and repeats the process until a valid factor is found. In the extremely unlikely case that Pollard's rho fails to find a valid solution, trial division is used as a backup, making this implementation 100% deterministic.

Trial division is always used for numbers below $2^{25}$, as the overhead of Pollard's rho makes it inefficient for smaller factors. See [performance](#performance-1) for more information.

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any number within the 64-bit integer range can be factorized in less than 21'000 nanoseconds, or 21 microseconds, with the majority of smaller numbers being a lot faster than that. The slowest cases occur near the upper end of the 64-bit range, where the factorization takes around 21 microseconds. The benchmark results can be found in [Documentation/Prime decomposition data.csv](./Documentation/Prime%20decomposition%20data.csv).

It also shows that trial division is faster for numbers up to about $2^{25}$, or roughly 30 million. To take advantage of that speed, a hybrid approach is used. Any factor below around $2^{25}$ is handled with trial division, while larger factors are found using Pollard’s rho.

<p align="center">
    <a href="./Documentation/Prime decomposition chart.svg">
        <img src="./Documentation/Prime decomposition chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
std::vector<q64bp::PrimeFactor> q64bp::prime_decomposition(q64bp::ui64 number);
```

#### Example

```c++
#include "Quick64BitPrimes/prime_decomposition.hpp"
#include <iostream>

int main() {

    int number = 123;

    auto prime_factors = q64bp::prime_decomposition(number);

    for (const auto& prime_factor : prime_factors) {

        std::cout << prime_factor.base << "^" << prime_factor.exponent << std::endl;

    }

}
```

## Tonelli-Shanks algorithm

To find the square root of a number modulo a prime, that is $r^{2} \equiv n \pmod p$, the Tonelli-Shanks algorithm can be used.

Tonelli-Shanks is a deterministic algorithm that finds the square root of a number modulo a prime by decomposing the prime's group structure and iteratively refining a candidate root using a non-square number until the square root is found. Its time complexity is typically $O(\log^{2}p)$.

This implementation uses fast paths wherever possible, often skipping the main iteration loop and significantly speeding up the algorithm. It generally returns both possible roots, $r$ and $p - r$, except in the trivial cases of $p = 2$ or $n = 0$, where only a single solution exists. In cases where no solution exists, nothing is returned. See [usage](#usage-2) for more information.

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any square root of a number modulo a prime within the 64-bit integer range can be found in less than 140 nanoseconds. The slowest cases occur near the upper end of the 64-bit range, where the algorithm takes around 130 nanoseconds. Most of the test duration is spent on the primality check. This is why the performance graph is almost identical to the one of the Miller-Rabin primality test. The benchmark results can be found in [Documentation/Tonelli-Shanks algorithm data.csv](./Documentation/Tonelli-Shanks%20algorithm%20data.csv).

<p align="center">
    <a href="./Documentation/Tonelli-Shanks algorithm chart.svg">
        <img src="./Documentation/Tonelli-Shanks algorithm chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
std::optional<std::pair<q64bp::ui64, std::optional<q64bp::ui64>>> q64bp::tonelli_shanks_algorithm(q64bp::ui64 number, q64bp::ui64 prime);
```

#### Example

```c++
#include "Quick64BitPrimes/tonelli_shanks_algorithm.hpp"
#include <iostream>

int main() {

    int number = 123;
    int prime = 137;

    auto modular_square_roots = q64bp::tonelli_shanks_algorithm(number, prime);

    if (modular_square_roots) {

        std::cout << modular_square_roots->first << std::endl;

        if (modular_square_roots->second) {

            std::cout << *modular_square_roots->second << std::endl;

        }

    }

}
```

## Fermat's sum of two squares theorem

To find the sum of two squares representation of a prime number, that is $x^{2} + y^{2} = p$, Fermat's sum of two squares theorem can be used.

This algorithm checks if the provided prime is in the form $p \equiv 1 \pmod 4$ and uses Tonelli-Shanks to get $r^{2} \equiv p-1 \pmod p$. A basic implementation of the Euclidean algorithm is used to find the two square roots of $r$, which are also the sum of two squares representation of the prime. Because this implementation leverages Tonelli-Shanks, its time complexity is approximately $O(\log^{2}p)$.

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any sum of two squares representation for any prime within the 64-bit integer range can be found in less than 140 nanoseconds. The slowest cases occur near the upper end of the 64-bit range, where the algorithm takes around 130 nanoseconds. Most of the test duration is spent on the primality check. This is why the performance graph is almost identical to the one of the Miller-Rabin primality test. The benchmark results can be found in [Documentation/Fermat's sum of two squares theorem data.csv](./Documentation/Fermat's%20sum%20of%20two%20squares%20theorem%20data.csv).

<p align="center">
    <a href="./Documentation/Fermat's sum of two squares theorem chart.svg">
        <img src="./Documentation/Fermat's sum of two squares theorem chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
std::optional<std::pair<q64bp::ui64, q64bp::ui64>> q64bp::fermat_sum_of_two_squares_theorem(q64bp::ui64 prime);
```

#### Example

```c++
#include "Quick64BitPrimes/fermat_sum_of_two_squares_theorem.hpp"
#include <iostream>

int main() {

    int prime = 137;

    auto square_roots = q64bp::fermat_sum_of_two_squares_theorem(prime);

    if (square_roots) {

        std::cout << square_roots->first << std::endl;
        std::cout << square_roots->second << std::endl;

    }

}
```

# Installation

## Requirements

This repository was developed and tested on Debian and requires the following packages to be installed:

- git
- build-essential

Install packages:
```sh
sudo apt install git build-essential
```

## Vendoring

The code in this repository is designed to be included as a library within an existing C++ project through vendoring. Instead of using a precompiled library, all header and source files are directly included and compiled within the project itself. A project including this library could look like this:

```
Project
├── Makefile
├── Source
│   └── main.cpp
└── Libraries
    └── Quick64BitPrimes
        ├── Include
        │   └── Quick64BitPrimes
        │       ├── fermat_sum_of_two_squares_theorem.hpp
        │       ├── miller_rabin_primality_test.hpp
        │       ├── modular_arithmetic.hpp
        │       ├── prime_decomposition.hpp
        │       ├── Quick64BitPrimes.hpp
        │       ├── tonelli_shanks_algorithm.hpp
        │       └── types.hpp
        └── Source
            ├── fermat_sum_of_two_squares_theorem.cpp
            ├── miller_rabin_primality_test.cpp
            ├── modular_arithmetic.cpp
            ├── prime_decomposition.cpp
            └── tonelli_shanks_algorithm.cpp
```

## Build

The code in this repository is designed to be included as a library within an existing C++ project through [vendoring](#vendoring). When compiled on its own, it produces an executable that accepts a single numeric input, determines whether the number is prime or returns its prime factors and reports the total execution time in nanoseconds.

Build the executable:
```sh
make
```

Run the executable:
```sh
./main
```

## Integer types and definitions

### `std::uint_fast64_t` & `q64bp::ui64`

This implementation uses the `std::uint_fast64_t` type wherever possible. This type guarantees at least 64 bits but may map to a larger and faster unsigned integer type if the platform provides one. The code aims to achieve optimal performance across different system architectures. The trade-off is that on systems where wider integer types are more efficient, the memory usage may increase slightly. The alias `q64bp::ui64` is defined in [Include/Quick64BitPrimes/types.hpp](./Include/Quick64BitPrimes/types.hpp) and can be changed from `std::uint_fast64_t` to `std::uint64_t` if a consistent memory footprint is preferred.

### `__uint128_t` & `q64bp::ui128`

The `__uint128_t` type is an unsigned 128-bit integer type available in GCC and Clang. It is used in parts of the modular arithmetic to safely handle large values without causing overflows. On compilers that do not support `__uint128_t`, a fallback method is used for modular multiplication. This fallback relies on repeated addition to avoid overflows and is typically around 10 times slower than using the 128-bit integer type. Therefore, availability of the `__uint128_t` type is critical for achieving optimal performance. The alias `q64bp::ui128` is defined in [Include/Quick64BitPrimes/types.hpp](./Include/Quick64BitPrimes/types.hpp).