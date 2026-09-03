// Compares modularfft::poly_multiply (NTT, exact modular arithmetic) against
// how convolution is usually done in practice: FFT over std::complex<double>.
// Two engines are benchmarked: a vendored dependency-free complex FFT
// (always built) and FFTW3 (built only if fftw3.h/libfftw3 are found at
// configure time; its column prints '-' otherwise).
//
// Prints three tables:
//   1. Speed, same operand sizes/values as bench_multiply.cpp.
//   2. Accuracy: small-digit coefficients (the realistic bignum-multiply
//      case) where modularfft's result is also the exact integer answer, so
//      it can serve as ground truth for the max rounding error of the
//      complex-FFT engines.
//   3. Scaling limits: since coefficients here are already `long long`, not
//      just digits, this sweeps coefficient magnitude against operand size
//      to find the smallest size at which each complex-FFT engine's first
//      rounding error appears, for each magnitude. Uses the largest usable
//      preset from data/primes.txt so the ceiling is modularFFT's own
//      transform-size limit, not an arbitrarily small modulus.
//
// Hand-rolled where possible; median of several runs per cell.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include <modularfft.hpp>

#include "complex_fft_ref.hpp"

#ifdef MODULARFFT_HAVE_FFTW3
#include <fftw3.h>
#endif

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

#ifdef MODULARFFT_HAVE_FFTW3
// Full plan-execute-destroy round trip per call, mirroring the vendored
// reference: each call is self-contained, no state reused across calls.
std::vector<std::int64_t> fftw_convolve(const std::vector<std::int64_t> &a,
                                        const std::vector<std::int64_t> &b) {
    if (a.empty() || b.empty()) {
        return {};
    }
    const std::size_t result_size = a.size() + b.size() - 1;
    const std::size_t n = complex_fft_ref::next_pow2(result_size);

    fftw_complex *fa = fftw_alloc_complex(n);
    fftw_complex *fb = fftw_alloc_complex(n);
    for (std::size_t i = 0; i < n; ++i) {
        fa[i][0] = i < a.size() ? static_cast<double>(a[i]) : 0.0;
        fa[i][1] = 0.0;
        fb[i][0] = i < b.size() ? static_cast<double>(b[i]) : 0.0;
        fb[i][1] = 0.0;
    }

    fftw_plan forward_a = fftw_plan_dft_1d(static_cast<int>(n), fa, fa, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan forward_b = fftw_plan_dft_1d(static_cast<int>(n), fb, fb, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_a);
    fftw_execute(forward_b);
    fftw_destroy_plan(forward_a);
    fftw_destroy_plan(forward_b);

    for (std::size_t i = 0; i < n; ++i) {
        const double re = fa[i][0] * fb[i][0] - fa[i][1] * fb[i][1];
        const double im = fa[i][0] * fb[i][1] + fa[i][1] * fb[i][0];
        fa[i][0] = re;
        fa[i][1] = im;
    }

    fftw_plan backward = fftw_plan_dft_1d(static_cast<int>(n), fa, fa, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    std::vector<std::int64_t> result(result_size);
    for (std::size_t i = 0; i < result_size; ++i) {
        result[i] = static_cast<std::int64_t>(std::llround(fa[i][0] / static_cast<double>(n)));
    }

    fftw_free(fa);
    fftw_free(fb);
    return result;
}
#endif

std::int64_t max_abs_error(const std::vector<modularfft::value_t> &exact,
                           const std::vector<std::int64_t> &approx) {
    std::int64_t worst = 0;
    for (std::size_t i = 0; i < exact.size(); ++i) {
        const std::int64_t diff = std::llabs(static_cast<std::int64_t>(exact[i]) - approx[i]);
        worst = std::max(worst, diff);
    }
    return worst;
}

} // namespace

int main() {
    constexpr modularfft::value_t kSizes[] = {256, 1024, 4096, 16384, 65536, 262144, 1048576};

    std::mt19937_64 rng(20260826);
    std::int64_t sink = 0; // keeps results from being optimised away

    std::cout << std::fixed << std::setprecision(3);

    std::cout << "-- speed: ntt vs complex-fft engines --\n";
    std::cout << "      n  ntt ms/op  ref ms/op  ref speedup  fftw ms/op  fftw speedup\n";
    for (modularfft::value_t n : kSizes) {
        const auto a =
            random_poly(rng, modularfft::preset_large.modulus, static_cast<std::size_t>(n));
        const auto b =
            random_poly(rng, modularfft::preset_large.modulus, static_cast<std::size_t>(n));

        const int repeats = n <= 16384 ? 7 : 5;
        const double ntt_ms =
            median_ms(repeats, [&] { sink += modularfft::poly_multiply(a, b).front(); });
        const double ref_ms =
            median_ms(repeats, [&] { sink += complex_fft_ref::convolve(a, b).front(); });

        std::cout << std::setw(7) << n << std::setw(11) << ntt_ms << std::setw(11) << ref_ms
                  << std::setw(12) << ref_ms / ntt_ms << 'x';

#ifdef MODULARFFT_HAVE_FFTW3
        const double fftw_ms = median_ms(repeats, [&] { sink += fftw_convolve(a, b).front(); });
        std::cout << std::setw(12) << fftw_ms << std::setw(13) << fftw_ms / ntt_ms << 'x';
#else
        std::cout << std::setw(12) << '-' << std::setw(14) << '-';
#endif
        std::cout << '\n';
    }

    // A wider coefficient range than the speed table's, and one size beyond
    // it: [0, 9] digits stay error-free out to several million points on
    // this machine (double precision has more headroom than intuition
    // suggests), so this range/size is chosen specifically to cross the
    // point where double-precision complex FFT starts mis-rounding —
    // demonstrating the actual failure mode NTT's exact arithmetic avoids.
    constexpr modularfft::value_t kAccuracySizes[] = {256,   1024,    4096,    16384,
                                                       65536, 262144, 1048576, 4194304};
    std::cout << "\n-- accuracy: max |error| vs exact NTT result, coefficients [0, 999] --\n";
    std::cout << "      n  ref max err  fftw max err\n";
    for (modularfft::value_t n : kAccuracySizes) {
        const auto a = random_poly(rng, 1000, static_cast<std::size_t>(n));
        const auto b = random_poly(rng, 1000, static_cast<std::size_t>(n));

        // Sums stay well under preset_large's modulus for these sizes and
        // coefficient range, so this result is also the exact integer answer.
        const auto exact = modularfft::poly_multiply(a, b, modularfft::preset_large);

        std::vector<std::int64_t> a64(a.begin(), a.end());
        std::vector<std::int64_t> b64(b.begin(), b.end());
        const auto ref_approx = complex_fft_ref::convolve(a64, b64);
        const std::int64_t ref_err = max_abs_error(exact, ref_approx);

        std::cout << std::setw(7) << n << std::setw(13) << ref_err;
#ifdef MODULARFFT_HAVE_FFTW3
        const auto fftw_approx = fftw_convolve(a64, b64);
        const std::int64_t fftw_err = max_abs_error(exact, fftw_approx);
        std::cout << std::setw(14) << fftw_err;
#else
        std::cout << std::setw(14) << '-';
#endif
        std::cout << '\n';
    }

    // Largest usable preset in data/primes.txt (modulus < 2^63, order_log2 =
    // 57 so transform-size is never the limiting factor here): coefficients
    // this library is actually built for are `long long`, not digits, so
    // this finds where each complex-FFT engine's first rounding error shows
    // up as magnitude climbs toward that range.
    const auto big_params = modularfft::make_params(4179340454199820289LL, 21, 57);
    constexpr modularfft::value_t kMagnitudes[] = {10,      100,       1000,      10000,
                                                    100000,  1000000,   10000000,
                                                    100000000, 1000000000};
    constexpr std::size_t kScalingSizes[] = {4,     8,     16,     32,     64,      128,
                                              256,   512,   1024,   2048,   4096,    8192,
                                              16384, 32768, 65536,  131072, 262144,  524288,
                                              1048576, 2097152, 4194304};

    std::cout << "\n-- scaling limits: first n with a rounding error, by coefficient "
                 "magnitude --\n";
    std::cout << "  magnitude   ref first-error n  fftw first-error n\n";
    for (modularfft::value_t magnitude : kMagnitudes) {
        std::size_t ref_first_error = 0;
        [[maybe_unused]] std::size_t fftw_first_error = 0;
        bool ref_done = false;
        bool fftw_done = false;
#ifndef MODULARFFT_HAVE_FFTW3
        fftw_done = true; // no column to fill in
#endif
        for (std::size_t n : kScalingSizes) {
            if (ref_done && fftw_done) {
                break;
            }
            // Stop once operands this large at this size would wrap
            // big_params' modulus themselves, since that would invalidate
            // the exact-result assumption used as ground truth below.
            const __int128 bound =
                static_cast<__int128>(n) * (magnitude - 1) * (magnitude - 1);
            if (bound >= static_cast<__int128>(big_params.modulus)) {
                break;
            }

            const auto a = random_poly(rng, magnitude, n);
            const auto b = random_poly(rng, magnitude, n);
            const auto exact = modularfft::poly_multiply(a, b, big_params);

            std::vector<std::int64_t> a64(a.begin(), a.end());
            std::vector<std::int64_t> b64(b.begin(), b.end());

            if (!ref_done) {
                const auto ref_approx = complex_fft_ref::convolve(a64, b64);
                if (max_abs_error(exact, ref_approx) > 0) {
                    ref_first_error = n;
                    ref_done = true;
                }
            }
#ifdef MODULARFFT_HAVE_FFTW3
            if (!fftw_done) {
                const auto fftw_approx = fftw_convolve(a64, b64);
                if (max_abs_error(exact, fftw_approx) > 0) {
                    fftw_first_error = n;
                    fftw_done = true;
                }
            }
#endif
        }

        std::cout << std::setw(11) << magnitude << std::setw(20)
                  << (ref_first_error ? std::to_string(ref_first_error) : "-") << std::setw(21);
#ifdef MODULARFFT_HAVE_FFTW3
        std::cout << (fftw_first_error ? std::to_string(fftw_first_error) : "-");
#else
        std::cout << '-';
#endif
        std::cout << '\n';
    }

    std::cout << "\n(checksum " << sink << ", ignore)\n";
    return 0;
}
