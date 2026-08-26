# modularFFT

Header-only C++17 library for the **modular fast Fourier transform** (a number-theoretic
transform over Z/p) and **O(n log n) polynomial multiplication** on top of it, based on
[Faster Integer Multiplication](https://www.csa.iisc.ac.in/~chandan/research/intMult.pdf)
by Daniel J. Bernstein.

All arithmetic is done modulo a prime `p` of the form `1 + k·2^i` with a known `2^i`-th root
of unity, using 128-bit intermediates for the products — so coefficients never overflow.

## Requirements

- GCC or Clang on a 64-bit platform (`__int128` support; a `#error` fires elsewhere)
- C++17 (only needed to *build* something that includes the header)

## Quick start

Copy [`include/modularfft.hpp`](include/modularfft.hpp) into your project, or use CMake:

```sh
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

```cpp
#include <modularfft.hpp>

modularfft::poly_multiply({1, 2, 3}, {4, 5});
// -> {4, 13, 22, 15}, i.e. (1 + 2x + 3x^2)(4 + 5x) mod p
```

A runnable version lives in [examples/multiply.cpp](examples/multiply.cpp).

## API

Everything lives in the `modularfft` namespace; values are `std::int64_t` in `[0, modulus)`.

| Function | Description |
| --- | --- |
| `value_t mul_mod(a, b, mod)` | `(a·b) % mod` without overflow. |
| `value_t mod_pow(base, exp, mod)` | Modular exponentiation by squaring. |
| `value_t mod_inverse(a, prime)` | Modular inverse via Fermat's little theorem. |
| `Params make_params(modulus, root, order_log2)` | Validates a parameter set (see below) and returns it. Throws `std::invalid_argument` if `root` does not have multiplicative order exactly `2^order_log2`. |
| `std::vector<value_t> ntt(data, params)` | Forward transform: coefficients → evaluations at the powers of `root`. Size must be a power of two ≤ `2^order_log2`; empty input yields empty output. |
| `std::vector<value_t> intt(data, params)` | Inverse transform; undoes `ntt` exactly. |
| `std::vector<value_t> poly_multiply(a, b, params = default_params)` | Polynomial product, lowest degree first. Result has exactly `a.size() + b.size() - 1` coefficients; an empty operand yields an empty result. Throws `std::invalid_argument` when the operands need a transform larger than `2^order_log2`. |

Built-in parameter sets:

| Name | modulus | root | order_log2 | max transform |
| --- | --- | --- | --- | --- |
| `preset_small` | 786433 | 8 | 17 | 2^17 points |
| `preset_medium` | 2013265921 | 137 | 27 | 2^27 points |
| `preset_large` = `default_params` | 7881299347898369 | 9 | 49 | 2^49 points |

Note that results are always **modulo the chosen prime** — pick a preset whose prime is
larger than any coefficient you expect (for raw integer products of digits < 10, e.g.
`n · 9²` must stay below the prime).

## Choosing your own parameters

`data/primes.txt` lists ready-made parameter sets, one per line:

```
Prime i j {root, 1/root} 1/2^i
```

where `Prime = 1 + j·2^i` is prime, `root` has multiplicative order `2^i` modulo `Prime`,
and the last two columns are derived values you do **not** need — `make_params(prime,
root, i)` computes them itself. Example: line `7881299347898369 49 14 {9, …} …` becomes

```cpp
auto params = modularfft::make_params(7881299347898369, 9, 49);
```

⚠️ Only rows whose `Prime` fits in a signed 64-bit integer (`Prime < 2^63`) are usable;
rows near the end of the file exceed that.

## Tests and benchmarks

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure   # unit tests
./build/bench_multiply                       # timing table
```

Tests check known products against hand-computed answers, random polynomials against a
schoolbook O(n²) reference, `intt(ntt(x)) == x` round-trips, and the error paths. The
benchmark compares NTT multiplication with schoolbook multiplication across operand sizes.
Each can be disabled with `MODULARFFT_BUILD_{EXAMPLES,TESTS,BENCHMARKS}=OFF`.

## License

[MIT](LICENSE)
