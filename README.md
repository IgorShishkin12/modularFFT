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
./build/bench_multiply                       # NTT vs schoolbook timing table
./build/bench_complex_fft                    # NTT vs complex-double FFT timing + accuracy
```

Tests check known products against hand-computed answers, random polynomials against a
schoolbook O(n²) reference, `intt(ntt(x)) == x` round-trips, and the error paths.
`bench_multiply` compares NTT multiplication with schoolbook multiplication across operand
sizes. `bench_complex_fft` compares it against how convolution is usually done in
practice — FFT over `std::complex<double>` — using a vendored dependency-free FFT (always
built) and, if `libfftw3-dev` is installed (`fftw3.h` and `libfftw3` found at configure
time; e.g. `apt install libfftw3-dev` / `brew install fftw`), FFTW3 too. Alongside speed it
reports max rounding error against modularFFT's exact result, since double-precision FFT
convolution eventually mis-rounds as operand size and coefficient magnitude grow, while the
NTT never does.

Each can be disabled with `MODULARFFT_BUILD_{EXAMPLES,TESTS,BENCHMARKS}=OFF`.

### Example results

From a dev machine (Release build, single-threaded throughout — numbers will vary on yours):

```
$ ./build/bench_multiply
      n  ntt ms/op  ns/coef  naive ms/op  speedup
    256      1.051 4106.820        1.023    0.973x
   1024      3.164 3090.325        5.001    1.580x
   4096      7.941 1938.786       78.798    9.923x
  16384     36.738 2242.313            -         -
  65536    172.248 2628.303            -         -
 262144    719.143 2743.313            -         -
1048576   3155.127 3008.963            -         -

$ ./build/bench_complex_fft
-- speed: ntt vs complex-fft engines --
      n  ntt ms/op  ref ms/op  ref speedup  fftw ms/op  fftw speedup
    256      0.494      0.075       0.151x       0.053        0.107x
   1024      2.191      0.293       0.134x       0.081        0.037x
   4096      7.989      1.455       0.182x       0.419        0.053x
  16384     35.753      6.779       0.190x       1.724        0.048x
  65536    159.153     32.130       0.202x       9.862        0.062x
 262144    689.319    144.156       0.209x      28.473        0.041x
1048576   3061.535    630.050       0.206x     169.970        0.056x

-- accuracy: max |error| vs exact NTT result, coefficients [0, 999] --
      n  ref max err  fftw max err
    256            0             0
   1024            0             0
   4096            0             0
  16384            0             0
  65536            0             0
 262144            0             0
1048576            4             0
4194304           96             0

-- scaling limits: first n with a rounding error, by coefficient magnitude --
  magnitude   ref first-error n  fftw first-error n
         10                   -                    -
        100             4194304                    -
       1000              524288                    -
      10000               32768                    -
     100000                8192               524288
    1000000                 512                 8192
   10000000                  32                  128
  100000000                   8                    4
 1000000000                   4                    4
```

The third table is the direct answer to "where does this stop scaling": coefficients here
are `long long`, not digits, and modularFFT's own ceiling comes from its modulus, not from
transform size — so the sweep uses the largest prime in `data/primes.txt` that still fits a
signed 64-bit integer (`4179340454199820289`, supporting transforms up to 2^57 points) and
finds, for each coefficient magnitude, the smallest operand size at which each complex-FFT
engine first mis-rounds a coefficient. Two things fall out of it:

- The crossover point collapses fast as magnitude grows. At magnitude 10 (single-digit) even
  4M-point transforms stayed exact; by magnitude 10⁹ — a plausible base-2^30-ish limb size
  for bignum arithmetic, and still far below what a 64-bit coefficient can hold — both
  engines fail at n = 4, the smallest size tested.
- FFTW's numerics are meaningfully better than the vendored radix-2 FFT (it holds out to a
  larger n at every magnitude where both eventually fail), but the trend is the same: past
  some scale, floating-point convolution is *always* wrong somewhere, regardless of
  implementation quality. NTT has no such ceiling — it's either usable at a given size (fits
  under the chosen prime) or it throws, never silently wrong.

Takeaways: this NTT implementation is currently *slower* than both complex-FFT engines at
every size tested — FFTW in particular is 15-30x faster single-threaded, so raw throughput
is not (yet) the pitch. What NTT buys you is exactness: the vendored complex FFT already
mis-rounds three-digit coefficients at n ≥ 2^20, while FFTW's more careful numerics held up
through n = 4194304 there. That gap collapses fast with coefficient magnitude, as the
scaling-limits table above shows.

## License

[MIT](LICENSE)
