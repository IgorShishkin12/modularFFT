#pragma once

// Minimal vendored radix-2 complex FFT, standing in for "how popular
// libraries do integer convolution via FFT": transform to the complex
// domain in double precision, pointwise multiply, transform back, round to
// the nearest integer. No external dependencies, used as the always-built
// baseline in bench_complex_fft.cpp.

#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace complex_fft_ref {

using Complex = std::complex<double>;

inline std::size_t next_pow2(std::size_t n) {
    std::size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

// In-place iterative Cooley-Tukey; `invert` selects forward/inverse.
inline void fft(std::vector<Complex> &a, bool invert) {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = (invert ? 1.0 : -1.0) * 2.0 * M_PI / static_cast<double>(len);
        const Complex wlen(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; ++k) {
                const Complex u = a[i + k];
                const Complex v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (Complex &x : a) {
            x /= static_cast<double>(n);
        }
    }
}

// Convolves two integer sequences via complex FFT, rounding each result
// coefficient to the nearest integer. Result has a.size() + b.size() - 1
// entries (empty if either input is empty).
inline std::vector<std::int64_t> convolve(const std::vector<std::int64_t> &a,
                                          const std::vector<std::int64_t> &b) {
    if (a.empty() || b.empty()) {
        return {};
    }
    const std::size_t result_size = a.size() + b.size() - 1;
    const std::size_t n = next_pow2(result_size);

    std::vector<Complex> fa(n, Complex(0.0, 0.0));
    std::vector<Complex> fb(n, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        fa[i] = Complex(static_cast<double>(a[i]), 0.0);
    }
    for (std::size_t i = 0; i < b.size(); ++i) {
        fb[i] = Complex(static_cast<double>(b[i]), 0.0);
    }

    fft(fa, false);
    fft(fb, false);
    for (std::size_t i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }
    fft(fa, true);

    std::vector<std::int64_t> result(result_size);
    for (std::size_t i = 0; i < result_size; ++i) {
        result[i] = static_cast<std::int64_t>(std::llround(fa[i].real()));
    }
    return result;
}

} // namespace complex_fft_ref
