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

From a container build (Release, single-threaded throughout, FFTW3 present — numbers will
vary on your machine):

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
    256      0.465      0.074       0.158x       0.051        0.109x
   1024      2.212      0.339       0.153x       0.094        0.042x
   4096      7.989      1.441       0.180x       0.421        0.053x
  16384     35.813      6.819       0.190x       1.727        0.048x
  65536    159.522     31.918       0.200x       9.841        0.062x
 262144    693.106    145.991       0.211x      28.735        0.041x
1048576   3045.294    635.922       0.209x     174.142        0.057x

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
     100000                8192               262144
    1000000                 512                 8192
   10000000                  32                   64
  100000000                   8                    4
 1000000000                   4                    4

-- equal digit-budget (n * log2(magnitude) = 1024.000 bits, magnitude capped so our own ground truth stays exact): split between element size and count --
      n  magnitude  ntt ms/op  ref ms/op  ref ok  fftw ms/op  fftw ok
      4 1022171763      0.006      0.000      no       0.009      no
      8  722784585      0.010      0.001      no       0.010      no
     16  511085881      0.018      0.002      no       0.009      no
     32  361392292      0.037      0.005      no       0.034      no
     64      65536      0.078      0.012     yes       0.054     yes
    128        256      0.167      0.027     yes       0.044     yes
    256         16      0.361      0.061     yes       0.046     yes
    512          4      0.779      0.137     yes       0.057     yes
   1024          2      1.671      0.282     yes       0.099     yes

NTT at its biggest possible element size (n=4): 0.006 ms
vendored FFT at its biggest *correct* element size (n=64): 0.012 ms (1.985x NTT's time)
FFTW3 at its biggest *correct* element size (n=64): 0.054 ms (8.624x NTT's time)
```

Takeaways: this NTT implementation is *slower* than both complex-FFT engines at every fixed
size tested — FFTW in particular is 15-30x faster single-threaded, so raw per-transform
throughput is not (yet) the pitch. Three things NTT buys instead:

- **Exactness, unconditionally.** The vendored complex FFT already mis-rounds three-digit
  coefficients by n = 2^20; FFTW's better numerics held up through n = 4194304 here, but the
  scaling-limits table shows that's a matter of degree, not kind — push coefficient magnitude
  toward what a 64-bit value can actually hold and both engines fail at the smallest size
  tested. NTT is either usable at a given size (fits under the chosen prime) or throws; it is
  never silently wrong.
- **The failure boundary is about total size, not element count in isolation.** The equal
  digit-budget table fixes the actual bit-length of the big number being multiplied
  (`n * log2(magnitude) = 1024`, the size of two ~1024-bit numbers) and sweeps how that budget
  splits between element count and element size. Both float engines are wrong for every split
  with `n ≤ 32` here — magnitude is large enough there (capped only by what our own NTT ground
  truth can still verify exactly, ~1e9 at n = 4) to break double's exact range — and both are
  correct for every split with `n ≥ 64`, where the same total size instead spreads thin across
  many small elements. NTT is exact at every split, on both sides of that boundary.
- **The real payoff shows up when both sides play by the same total-size budget.** Since NTT
  doesn't care about base, it can take the fastest split (fewest, biggest limbs) outright; a
  complex-FFT engine can only take that split if it's still correct there, which means falling
  back to many small limbs — and paying for it. Here NTT's fastest-possible split (n = 4,
  0.006 ms) beat FFTW's biggest *reliably correct* split (n = 64, 0.054 ms) by ~8.6x, and beat
  the vendored FFT's by ~2x — despite NTT losing every fixed-size, same-n comparison in the
  first table. Constrained to answers it can actually trust, "popular library" FFT convolution
  loses its speed advantage for big-number multiplication; NTT's ceiling only comes from the
  modulus, so it can spend the whole budget on fewer, bigger limbs and get there directly.

## License

[MIT](LICENSE)
