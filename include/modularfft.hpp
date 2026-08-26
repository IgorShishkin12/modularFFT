#pragma once

// Header-only modular fast Fourier transform (NTT) over Z/p, with O(n log n)
// polynomial multiplication on top. Based on "Faster Integer Multiplication"
// by Dan Bernstein, https://www.csa.iisc.ac.in/~chandan/research/intMult.pdf
//
// Requires a 64-bit compiler with __int128 support (GCC or Clang).

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#if !defined(__SIZEOF_INT128__)
#error "modularfft needs __int128 support (compile with GCC or Clang on 64-bit)"
#endif

namespace modularfft {

// Coefficients and transformed values; always interpreted modulo `modulus`.
using value_t = std::int64_t;

// (a * b) % mod via a 128-bit intermediate, so no overflow for mod < 2^63.
[[nodiscard]] inline value_t mul_mod(value_t a, value_t b, value_t mod) {
    return static_cast<value_t>((static_cast<__int128>(a) * static_cast<__int128>(b)) %
                                static_cast<__int128>(mod));
}

// base^exp % mod by exponentiation by squaring.
[[nodiscard]] inline value_t mod_pow(value_t base, value_t exp, value_t mod) {
    value_t result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// Modular inverse of a modulo a prime, via Fermat's little theorem.
[[nodiscard]] inline value_t mod_inverse(value_t a, value_t prime) {
    return mod_pow(a, prime - 2, prime);
}

// Transform parameters over Z/modulus.
//
// `modulus` must be a prime of the form 1 + k*2^order_log2, and `root` an
// element of multiplicative order exactly 2^order_log2 (a 2^order_log2-th
// root of unity). Transforms of any power-of-two size up to 2^order_log2 are
// supported; the inverse root and normalisation factor are derived
// internally. See data/primes.txt for ready-made parameter sets.
struct Params {
    value_t modulus;
    value_t root;
    std::size_t order_log2;
};

// Validates and packages transform parameters. Throws std::invalid_argument
// if they cannot describe a 2^order_log2-th root of unity. Note this checks
// the order of `root`; proving that `modulus` is actually prime is left to
// the caller.
[[nodiscard]] inline Params make_params(value_t modulus, value_t root, std::size_t order_log2) {
    if (modulus < 3 || order_log2 == 0 || order_log2 > 62) {
        throw std::invalid_argument("need a prime modulus >= 3 and 0 < order_log2 <= 62");
    }
    if (root <= 1 || root >= modulus) {
        throw std::invalid_argument("root must lie strictly between 1 and modulus");
    }
    const value_t order = static_cast<value_t>(1) << order_log2;
    if (mod_pow(root, order, modulus) != 1 || mod_pow(root, order / 2, modulus) == 1) {
        throw std::invalid_argument("root does not have multiplicative order 2^order_log2");
    }
    return Params{modulus, root, order_log2};
}

// Ready-made parameter sets, taken from data/primes.txt.
inline constexpr Params preset_small{786433, 8, 17};           // transforms up to 2^17 points
inline constexpr Params preset_medium{2013265921, 137, 27};    // up to 2^27 points
inline constexpr Params preset_large{7881299347898369, 9, 49}; // up to 2^49 points
inline constexpr Params default_params = preset_large;

namespace detail {

inline std::size_t next_power_of_two(std::size_t x) {
    std::size_t power = 1;
    while (power < x) {
        power <<= 1;
    }
    return power;
}

inline std::size_t log2_of_power_of_two(std::size_t power_of_two) {
    std::size_t log = 0;
    while ((std::size_t{1} << log) < power_of_two) {
        ++log;
    }
    return log;
}

// Checks that n is a usable transform length for `params` (non-zero power of
// two, at most 2^order_log2) and squares `root` down from order 2^order_log2
// to order exactly n.
inline value_t scaled_root(const Params &params, std::size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw std::invalid_argument("transform size must be a non-zero power of two");
    }
    if (params.order_log2 > 62 || n > (std::size_t{1} << params.order_log2)) {
        throw std::invalid_argument("transform size exceeds 2^order_log2 for these params");
    }
    const std::size_t target_order = log2_of_power_of_two(n);
    value_t root = params.root;
    for (std::size_t order = params.order_log2; order > target_order; --order) {
        root = mul_mod(root, root, params.modulus); // halve the root's order
    }
    return root;
}

// Recursive Cooley-Tukey step: evaluates `data` (coefficients, lowest degree
// first) at the powers of `root`, where root must have order exactly
// data.size(). Sequences of length <= 1 are their own transform.
[[nodiscard]] inline std::vector<value_t> transform(std::vector<value_t> data, value_t root,
                                                    value_t modulus) {
    const std::size_t n = data.size();
    if (n <= 1) {
        return data;
    }

    // Transform the even- and odd-indexed halves with the squared root
    // (which has half the order).
    std::vector<value_t> even(n / 2), odd(n / 2);
    for (std::size_t i = 0; i < n; i += 2) {
        even[i / 2] = data[i];
        odd[i / 2] = data[i + 1];
    }
    const value_t root_squared = mul_mod(root, root, modulus);
    even = transform(std::move(even), root_squared, modulus);
    odd = transform(std::move(odd), root_squared, modulus);

    // Combine: X[k] = E[k mod n/2] + root^k * O[k mod n/2].
    std::vector<value_t> result(n);
    const std::size_t half = n / 2;
    value_t w = 1; // w = root^k
    for (std::size_t k = 0; k < n; ++k) {
        result[k] = (even[k % half] + mul_mod(w, odd[k % half], modulus)) % modulus;
        w = mul_mod(w, root, modulus);
    }
    return result;
}

} // namespace detail

// Forward transform: maps coefficients to evaluations at the powers of the
// parameter root. `data.size()` must be a power of two not exceeding
// 2^order_log2 (an empty input yields an empty output), and coefficients
// must already lie in [0, modulus).
[[nodiscard]] inline std::vector<value_t> ntt(std::vector<value_t> data, const Params &params) {
    if (data.empty()) {
        return data;
    }
    const value_t root = detail::scaled_root(params, data.size());
    return detail::transform(std::move(data), root, params.modulus);
}

// Inverse transform: undoes ntt exactly (scaling by 1/n internally).
[[nodiscard]] inline std::vector<value_t> intt(std::vector<value_t> data, const Params &params) {
    if (data.empty()) {
        return data;
    }
    Params inverse = params;
    inverse.root = mod_inverse(params.root, params.modulus);
    data = ntt(std::move(data), inverse);
    const value_t scale = mod_inverse(static_cast<value_t>(data.size()), params.modulus);
    for (value_t &coefficient : data) {
        coefficient = mul_mod(coefficient, scale, params.modulus);
    }
    return data;
}

// Multiplies polynomials `a` and `b` (coefficients lowest degree first,
// values in [0, modulus)) in O(n log n). The result holds exactly
// a.size() + b.size() - 1 coefficients; an empty operand yields an empty
// product. Throws std::invalid_argument if the operands require a transform
// larger than 2^order_log2 — switch to a parameter set with a bigger order.
[[nodiscard]] inline std::vector<value_t> poly_multiply(const std::vector<value_t> &a,
                                                        const std::vector<value_t> &b,
                                                        const Params &params = default_params) {
    if (a.empty() || b.empty()) {
        return {};
    }

    const std::size_t product_size = a.size() + b.size() - 1;
    const std::size_t n = detail::next_power_of_two(product_size);

    std::vector<value_t> fa(a), fb(b);
    fa.resize(n, 0);
    fb.resize(n, 0);
    fa = ntt(std::move(fa), params);
    fb = ntt(std::move(fb), params);

    for (std::size_t i = 0; i < n; ++i) {
        fa[i] = mul_mod(fa[i], fb[i], params.modulus);
    }

    std::vector<value_t> product = intt(std::move(fa), params);
    product.resize(product_size); // strip the zero padding beyond the true degree
    return product;
}

} // namespace modularfft
