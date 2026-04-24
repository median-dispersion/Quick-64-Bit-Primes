# Quick 64-Bit Primes

This repository contains a C++ implementation of the [Miller-Rabin primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test) and [Pollard's rho prime decomposition / factorization algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm) using [Brent's cycle detection](https://en.wikipedia.org/wiki/Cycle_detection#Brent.27s_algorithm) method for the full range of unsigned 64-bit integers. It accepts a number in between $0$ and $2^{64}-1$ and, depending on the chosen method, either determines whether the number is prime or returns all of its prime factors. The runtime ranges from just a few nanoseconds up to about 20 microseconds in the worst case for a full 64-bit value.

## Primality test

To quickly check if a number is prime, the Miller-Rabin primality test can be used.

### Probabilistic algorithm, deterministic result

The Miller-Rabin primality test is a probabilistic algorithm that checks whether a number is prime by testing if it behaves like a prime under repeated modular exponentiation. If it fails any test, it’s guaranteed to be composite, but if it passes multiple rounds, it’s very likely prime. Its time complexity is $O(k\log^{3}n)$, where $n$ is the number being tested and $k$ is the number of rounds.

Probabilistic, in this context, means the test does not always produce a 100% correct result and may occasionally incorrectly classify a number as prime. However, this particular implementation is both deterministic and fully accurate. Restricting the input to the unsigned 64-bit integer range allows the use of a fixed set of bases that guarantees 100% correct primality testing within that range. It also simplifies the complexity to $O(\log^{3}n)$.

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any number within the 64-bit integer range can be checked in less than 100 nanoseconds. The slowest cases occur near the upper end of the 64-bit range, where the test takes around 100 nanoseconds. The benchmark results can be found in [Documentation/Primality test data.csv](./Documentation/Primality%20test%20data.csv).

<p align="center">
    <a href="./Documentation/Primality test chart.svg">
        <img src="./Documentation/Primality test chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
bool q64bp::millerRabinPrimalityTest(std::uint_fast64_t number);
```

#### Example

```c++
#include "Quick64BitPrimes.hpp"
#include <iostream>

int main() {

    int number = 123;

    if (q64bp::millerRabinPrimalityTest(number)) {

        std::cout << number << " is prime!" << std::endl;

    } else {

        std::cout << number << " is not prime." << std::endl;

    }

}
```

## Prime decomposition

To decompose a number into its prime factors, Brent's improved variant of Pollard's rho can be used.

### Probabilistic algorithm, deterministic result

Pollard's rho is a probabilistic algorithm that finds a nontrivial factor of a composite number $n$ by generating a pseudo-random sequence modulo $n$ and using cycle detection to reveal a hidden common divisor. Brent's variant of Pollard's rho is a cycle-detection improvement over Floyd's method that batches iterations to reduce the expensive greatest common divisor computations while still finding a nontrivial factor of a composite number. Its worst-case time complexity is $O(n^{1/4})$, and in practice it's often faster.

Probabilistic, in this context, means that the algorithm might fail to find a factor that is not one or the number itself. However, this particular implementation of Pollard’s rho behaves deterministically in practice. Whenever it fails, it simply retries with a new set of random values for the polynomial function and repeats the process, exhausting all possible combinations, until a valid factor is found, making it "deterministic".

### Performance

All tests were performed single-threaded on an AMD Ryzen™ 5 7600X.

The benchmark results show that almost any number within the 64-bit integer range can be factorized in less than 20'000 nanoseconds, or 20 microseconds, with the majority of smaller numbers being a lot faster than that. The slowest cases occur near the upper end of the 64-bit range, where the factorization takes around 20 microseconds. The benchmark results can be found in [Documentation/Prime decomposition data.csv](./Documentation/Prime%20decomposition%20data.csv).

It also shows that trial division is faster for numbers up to about $2^{25}$, or roughly 30 million. To take advantage of that speed, a hybrid approach is used. Any factor below around $2^{25}$ is handled with trial division, while larger factors are found using Pollard’s rho.

<p align="center">
    <a href="./Documentation/Prime decomposition chart.svg">
        <img src="./Documentation/Prime decomposition chart.svg" style="width: 800px;">
    </a>
</p>

### Usage

#### Function prototype

```c++
std::vector<q64bp::PrimeFactor> q64bp::primeDecomposition(std::uint_fast64_t number);
```

#### Example

```c++
#include "Quick64BitPrimes.hpp"
#include <vector>
#include <iostream>

int main() {

    int number = 123;

    std::vector<q64bp::PrimeFactor> primeFactors = q64bp::primeDecomposition(number);

    for (auto& primeFactor : primeFactors) {

        std::cout << primeFactor.base << "^" << primeFactor.exponent << std::endl;

    }

}
```

## Requirements

This repository was developed and tested on Debian and requires the following packages to be installed:

- git
- g++
- make

Install packages:
```bash
sudo apt install git build-essential
```

## Build

The code is designed to be included as a library within an existing C++ project and used in that context. When compiled on its own, it produces an executable that accepts a single numeric input, determines whether the number is prime or returns its prime factors, and reports the total execution time in nanoseconds.

Clone the repository:
```bash
git clone https://github.com/median-dispersion/Quick-64-Bit-Primes.git
```

Enter the repository:
```bash
cd Quick-64-Bit-Primes
```

Build the repository:
```bash
make
```

Run the executable:
```bash
./main
```