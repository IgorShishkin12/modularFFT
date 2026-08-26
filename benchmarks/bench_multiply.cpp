// Times poly_multiply across operand sizes and compares it to schoolbook
// multiplication where that is still affordable. Prints one row per size:
//
//   n  ntt ms/op  ns/coef  naive ms/op  speedup
//
// Hand-rolled (no dependencies); median of several runs per cell.

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <modularfft.hpp>

namespace {

using Clock = std::chrono::steady_clock;

template <typename Fn> double median_ms(int repeats, Fn &&fn) {
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repeats));
    for (int i = 0; i < repeats; ++i) {
        const Clock::time_point start = Clock::now();
        fn();
        const Clock::time_point stop = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(stop - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

std::vector<modularfft::value_t> random_poly(std::mt19937_64 &rng, modularfft::value_t modulus,
                                             std::size_t length) {
    std::uniform_int_distribution<modularfft::value_t> dist(0, modulus - 1);
    std::vector<modularfft::value_t> poly(length);
    for (modularfft::value_t &coefficient : poly) {
        coefficient = dist(rng);
    }
    return poly;
}

// Schoolbook O(n^2) reference for comparison.
modularfft::value_t naive_checksum(const std::vector<modularfft::value_t> &a,
                                   const std::vector<modularfft::value_t> &b,
                                   modularfft::value_t modulus) {
    modularfft::value_t sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            sum += modularfft::mul_mod(a[i], b[j], modulus);
        }
    }
    return sum % modulus;
}

} // namespace

int main() {
    constexpr modularfft::value_t kSizes[] = {256, 1024, 4096, 16384, 65536, 262144, 1048576};

    std::mt19937_64 rng(20260826);
    modularfft::value_t sink = 0; // keeps the products from being optimised away

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "      n  ntt ms/op  ns/coef  naive ms/op  speedup\n";
    for (modularfft::value_t n : kSizes) {
        const auto a =
            random_poly(rng, modularfft::preset_large.modulus, static_cast<std::size_t>(n));
        const auto b =
            random_poly(rng, modularfft::preset_large.modulus, static_cast<std::size_t>(n));

        const int repeats = n <= 16384 ? 7 : 5;
        const double ntt_ms =
            median_ms(repeats, [&] { sink += modularfft::poly_multiply(a, b).front(); });

        const bool naive_affordable = n <= 4096;
        double naive_ms = 0.0;
        if (naive_affordable) {
            naive_ms = median_ms(
                repeats, [&] { sink += naive_checksum(a, b, modularfft::preset_large.modulus); });
        }

        std::cout << std::setw(7) << n << std::setw(11) << ntt_ms << std::setw(9)
                  << ntt_ms * 1e6 / static_cast<double>(n);
        if (naive_affordable) {
            std::cout << std::setw(13) << naive_ms << std::setw(9) << naive_ms / ntt_ms << 'x';
        } else {
            std::cout << std::setw(13) << '-' << std::setw(10) << '-';
        }
        std::cout << '\n';
    }

    std::cout << "\n(checksum " << sink << ", ignore)\n";
    return 0;
}
