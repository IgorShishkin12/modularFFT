// Multiplies 1 + 2x + 3x^2 by itself, with the default parameters and with
// a custom parameter set.
#include <iostream>
#include <vector>

#include <modularfft.hpp>

void print(const std::vector<modularfft::value_t>& poly) {
    for (modularfft::value_t coefficient : poly)
        std::cout << coefficient << ' ';
    std::cout << '\n';
}

int main() {
    const std::vector<modularfft::value_t> a{1, 2, 3};  // 1 + 2*x + 3*x^2

    print(modularfft::poly_multiply(a, a));
    // -> 1 4 10 12 9  (1 + 4*x + 10*x^2 + 12*x^3 + 9*x^4)

    // The same product with different parameters (row 786433 of
    // data/primes.txt). Everything is taken modulo the chosen prime.
    const modularfft::Params params = modularfft::make_params(786433, 8, 17);
    print(modularfft::poly_multiply(a, a, params));

    return 0;
}
