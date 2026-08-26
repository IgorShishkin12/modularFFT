// Tests for modularfft: known answers, property tests against schoolbook
// multiplication, transform round-trips and error handling.
//
// Hand-rolled runner (no dependencies): prints each failure and exits
// non-zero if any check failed.

#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <modularfft.hpp>

namespace {

int checks_run = 0;
int checks_failed = 0;

void expect_equal(const std::vector<modularfft::value_t>& actual,
                  const std::vector<modularfft::value_t>& expected,
                  const std::string& what) {
    ++checks_run;
    if (actual != expected) {
        ++checks_failed;
        std::cerr << "FAIL: " << what << "\n  expected:";
        for (modularfft::value_t v : expected)
            std::cerr << ' ' << v;
        std::cerr << "\n  actual:  ";
        for (modularfft::value_t v : actual)
            std::cerr << ' ' << v;
        std::cerr << '\n';
    }
}

// Schoolbook O(n^2) reference, independent of the NTT code path.
std::vector<modularfft::value_t> naive_multiply(const std::vector<modularfft::value_t>& a,
                                                const std::vector<modularfft::value_t>& b,
                                                modularfft::value_t modulus) {
    if (a.empty() || b.empty())
        return {};
    std::vector<modularfft::value_t> product(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            product[i + j] = (product[i + j] + modularfft::mul_mod(a[i], b[j], modulus)) % modulus;
    return product;
}

std::vector<modularfft::value_t> random_poly(std::mt19937_64& rng,
                                             modularfft::value_t modulus,
                                             std::size_t length,
                                             bool include_boundaries) {
    std::uniform_int_distribution<modularfft::value_t> dist(0, modulus - 1);
    std::vector<modularfft::value_t> poly(length);
    for (modularfft::value_t& coefficient : poly)
        coefficient = dist(rng);
    if (include_boundaries && !poly.empty()) {
        poly.front() = modulus - 1;  // exercise wrap-around
        poly.back() = 0;
    }
    return poly;
}

void test_known_answers() {
    using modularfft::poly_multiply;
    using modularfft::value_t;

    const std::vector<value_t> cubic{1, 2, 3};  // 1 + 2*x + 3*x^2
    const std::vector<value_t> square{1, 4, 10, 12, 9};
    expect_equal(poly_multiply(cubic, cubic), square, "(1+2x+3x^2)^2 with default params");
    expect_equal(poly_multiply(cubic, cubic, modularfft::preset_small), square,
                 "(1+2x+3x^2)^2 with preset_small");

    expect_equal(poly_multiply({0, 0}, {5, 7}), {0, 0, 0}, "zero polynomial");
    expect_equal(poly_multiply({12}, {4}), {48}, "constants");
    expect_equal(poly_multiply({2}, {2, 3}), {4, 6}, "constant times linear");
    expect_equal(poly_multiply({}, {1, 2}), {}, "empty operand");

    // (p-1)^2 == 1 (mod p) for any of the presets.
    expect_equal(poly_multiply({modularfft::default_params.modulus - 1},
                               {modularfft::default_params.modulus - 1}),
                 {1}, "(p-1)^2 mod p");
    expect_equal(poly_multiply({modularfft::preset_small.modulus - 1},
                               {modularfft::preset_small.modulus - 1},
                               modularfft::preset_small),
                 {1}, "(p_small-1)^2 mod p_small");
}

void test_transform_structure() {
    using namespace modularfft;

    // The constant polynomial 1 evaluates to 1 at every power of the root.
    expect_equal(ntt({1, 0, 0, 0}, preset_small), {1, 1, 1, 1},
                 "ntt of delta_0 is all ones");

    // Round-trip: intt(ntt(x)) == x for a spread of sizes and both presets.
    std::mt19937_64 rng(20260826);
    for (std::size_t size = 1; size <= 4096; size *= 2) {
        for (const Params* params : {&preset_small, &preset_large}) {
            const std::vector<value_t> x = random_poly(rng, params->modulus, size, false);
            expect_equal(intt(ntt(x, *params), *params), x,
                         "intt(ntt(x)) round-trip, size " + std::to_string(size));
        }
    }

    expect_equal(ntt({}, preset_large), {}, "ntt of empty input");
    expect_equal(intt({}, preset_large), {}, "intt of empty input");
}

void test_property_against_naive() {
    using namespace modularfft;

    static constexpr std::size_t lengths[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32,
                                              33, 63, 64, 65, 100, 127, 128, 129, 200, 255,
                                              256, 257, 300};
    std::mt19937_64 rng(1234567);

    int iteration = 0;
    for (std::size_t length_a : lengths) {
        std::uniform_int_distribution<std::size_t> length_dist(1, 300);
        const std::size_t length_b = length_dist(rng);
        const bool boundaries = (iteration++ % 4 == 0);

        for (const Params* params : {&preset_small, &preset_large}) {
            const std::vector<value_t> a = random_poly(rng, params->modulus, length_a, boundaries);
            const std::vector<value_t> b = random_poly(rng, params->modulus, length_b, boundaries);
            expect_equal(poly_multiply(a, b, *params),
                         naive_multiply(a, b, params->modulus),
                         "random product vs schoolbook, sizes " + std::to_string(length_a) +
                             "+" + std::to_string(length_b));
        }
    }
}

template <typename Fn>
void expect_invalid_argument(Fn&& fn, const std::string& what) {
    ++checks_run;
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return;
    } catch (...) {
        ++checks_failed;
        std::cerr << "FAIL: " << what << " threw something other than std::invalid_argument\n";
        return;
    }
    ++checks_failed;
    std::cerr << "FAIL: " << what << " did not throw\n";
}

void test_error_handling() {
    using namespace modularfft;

    // preset_small supports transforms up to 2^17 points; this needs 2^18.
    const std::vector<value_t> big(100000, 1);
    expect_invalid_argument([&] { (void)poly_multiply(big, big, preset_small); },
                            "operands beyond capacity throw");

    expect_invalid_argument([] { (void)ntt({1, 2, 3}, preset_large); },
                            "non-power-of-two ntt size throws");

    expect_invalid_argument([] { (void)make_params(786433, 1, 17); }, "root == 1 rejected");
    expect_invalid_argument([] { (void)make_params(786433, 0, 17); }, "root == 0 rejected");
    expect_invalid_argument([] { (void)make_params(786433, 8, 16); },
                            "wrong-order root rejected");
    expect_invalid_argument([] { (void)make_params(15, 2, 3); },
                            "order check rejects non-prime modulus");
}

}  // namespace

int main() {
    test_known_answers();
    test_transform_structure();
    test_property_against_naive();
    test_error_handling();

    std::cout << checks_run << " checks, " << checks_failed << " failed\n";
    return checks_failed == 0 ? 0 : 1;
}
