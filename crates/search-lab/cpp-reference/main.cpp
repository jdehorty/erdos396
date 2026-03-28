// erdos396-speed — Minimal, highly-optimized witness search for Erdős #396
// Ported from third_party/erdos_problem_396.cpp.
// Build: make (see Makefile for flags)

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,unroll-loops")
#endif
#if defined(__GNUC__) && defined(__x86_64__) && !defined(__clang__)
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#endif

#include <algorithm>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

// ---------------------------------------------------------------------------
// Tunables
// ---------------------------------------------------------------------------
static constexpr uint64_t CHUNK_SIZE = 10'000'000;
static constexpr uint64_t BLOCK_SIZE = 32768;   // fits L1 cache
static constexpr uint64_t NUM_SEQ    = 2;

static unsigned g_threads = std::thread::hardware_concurrency();
static std::vector<uint64_t> g_primes;

// ---------------------------------------------------------------------------
// Sieve of Eratosthenes
// ---------------------------------------------------------------------------
static void sieve_primes(uint64_t limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (uint64_t p = 2; p * p <= limit; ++p)
        if (is_prime[p])
            for (uint64_t m = p * p; m <= limit; m += p)
                is_prime[m] = false;
    g_primes.reserve(limit / 10);
    for (uint64_t p = 2; p <= limit; ++p)
        if (is_prime[p])
            g_primes.push_back(p);
}

// ---------------------------------------------------------------------------
// Exact divisibility check  (Legendre formula, for final validation)
// ---------------------------------------------------------------------------
static bool exact_check(uint64_t n, uint64_t k) {
    uint64_t nu2_prod = std::popcount(n - k - 1) - std::popcount(n) + k + 1;
    if (static_cast<uint64_t>(std::popcount(n)) < nu2_prod) [[unlikely]]
        return false;

    for (uint64_t p : g_primes) {
        if (p == 2) continue;
        if (p > 2 * k) break;
        uint64_t nu_prod = 0, nu_comb = 0, pw = p;
        for (;;) {
            uint64_t vn = n / pw;
            nu_prod += vn - (n - k - 1) / pw;
            nu_comb += (2 * n) / pw - 2 * vn;
            if (pw > (2 * n) / p) break;
            pw *= p;
        }
        if (nu_prod > nu_comb) [[unlikely]]
            return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Template-specialized prime processor  (compiler turns / and % into
// multiply-shift for compile-time-known p)
// ---------------------------------------------------------------------------
template<uint64_t p>
static inline void process_prime(
    uint64_t start_c, uint64_t start_j, uint64_t W,
    uint8_t* __restrict__ vld, uint64_t* __restrict__ rem, uint64_t k)
{
    if constexpr (p > 0) {  // always true; keeps compiler happy
        if (p > 2 * k) {
            uint64_t c = start_c;
            for (uint64_t j = start_j; j < W; j += p, ++c) {
                if (!vld[j]) continue;
                bool ok = true;
                if (c * 2 < p) {
                    ok = false;
                } else if (c >= p) {
                    // Kummer carry check: v_p(C(2n,n)) > v_p(n)?
                    uint64_t nu = 0, t = c;
                    while (t > 0 && t % p == 0) { ++nu; t /= p; }
                    uint64_t carries = 0, carry = 0;
                    while (t > 0) {
                        uint64_t d = t % p;
                        if (d * 2 + carry >= p) {
                            ++carries; carry = 1;
                            if (carries > nu) break;
                        } else {
                            carry = 0;
                        }
                        t /= p;
                    }
                    if (carries <= nu) ok = false;
                }
                if (!ok) [[unlikely]] {
                    vld[j] = 0;
                } else {
                    uint64_t t = rem[j] / p;
                    if (t > 0) while (t % p == 0) t /= p;
                    rem[j] = t;
                }
            }
        } else {
            for (uint64_t j = start_j; j < W; j += p) {
                if (!vld[j]) continue;
                uint64_t t = rem[j] / p;
                if (t > 0) while (t % p == 0) t /= p;
                rem[j] = t;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Dynamic (non-specialized) prime processor  (primes > 199)
// ---------------------------------------------------------------------------
static inline void process_prime_dyn(
    uint64_t p, uint64_t start_c, uint64_t start_j, uint64_t W,
    uint8_t* __restrict__ vld, uint64_t* __restrict__ rem, uint64_t k)
{
    if (p > 2 * k) {
        uint64_t c = start_c;
        for (uint64_t j = start_j; j < W; j += p, ++c) {
            if (!vld[j]) continue;
            bool ok = true;
            if (c * 2 < p) {
                ok = false;
            } else if (c >= p) {
                uint64_t nu = 0, t = c;
                while (t > 0 && t % p == 0) { ++nu; t /= p; }
                uint64_t carries = 0, carry = 0;
                while (t > 0) {
                    uint64_t d = t % p;
                    if (d * 2 + carry >= p) {
                        ++carries; carry = 1;
                        if (carries > nu) break;
                    } else {
                        carry = 0;
                    }
                    t /= p;
                }
                if (carries <= nu) ok = false;
            }
            if (!ok) [[unlikely]] {
                vld[j] = 0;
            } else {
                uint64_t t = rem[j] / p;
                if (t > 0) while (t % p == 0) t /= p;
                rem[j] = t;
            }
        }
    } else {
        for (uint64_t j = start_j; j < W; j += p) {
            if (!vld[j]) continue;
            uint64_t t = rem[j] / p;
            if (t > 0) while (t % p == 0) t /= p;
            rem[j] = t;
        }
    }
}

// ---------------------------------------------------------------------------
// Compile-time dispatch macro
// ---------------------------------------------------------------------------
#define P(v) case v: process_prime<v>(sc, sj, W, vld, rm, k); break;

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------
static uint64_t solve(uint64_t k, uint64_t start_L) {
    const uint64_t cpb = g_threads * NUM_SEQ;
    uint64_t L_batch = start_L;

    std::atomic<uint64_t> best{UINT64_MAX};
    std::atomic<bool> found{false};

    auto worker = [&](std::atomic<int>& next_chunk) {
        std::vector<uint8_t>  valid_buf(BLOCK_SIZE);
        std::vector<uint64_t> rem_buf(BLOCK_SIZE);

        while (true) {
            int cid = next_chunk.fetch_add(1, std::memory_order_relaxed);
            if (cid >= static_cast<int>(cpb)) break;

            uint64_t L = L_batch + static_cast<uint64_t>(cid) * CHUNK_SIZE;
            uint64_t R = L + CHUNK_SIZE + k - 1;
            if (L > best.load(std::memory_order_relaxed)) continue;

            uint64_t max_p = std::max<uint64_t>(
                static_cast<uint64_t>(std::sqrt(2.0 * R)) + 1, 2 * k);
            uint64_t consec = 0;

            for (uint64_t bL = L; bL <= R; bL += BLOCK_SIZE) {
                uint64_t bR = std::min(R, bL + BLOCK_SIZE - 1);
                uint64_t W  = bR - bL + 1;
                auto* vld = valid_buf.data();
                auto* rm  = rem_buf.data();

                for (uint64_t j = 0; j < W; ++j) { rm[j] = bL + j; vld[j] = 1; }

                for (uint64_t p : g_primes) {
                    if (p > max_p) break;
                    uint64_t sc = (bL + p - 1) / p;
                    uint64_t sj = sc * p - bL;
                    if (sj >= W) continue;

                    if (p == 2) {
                        for (uint64_t j = sj; j < W; j += 2)
                            if (vld[j]) { uint64_t t = rm[j]; t >>= std::countr_zero(t); rm[j] = t; }
                        continue;
                    }

                    switch (p) {
                        P(3)   P(5)   P(7)   P(11)  P(13)  P(17)  P(19)  P(23)
                        P(29)  P(31)  P(37)  P(41)  P(43)  P(47)  P(53)  P(59)
                        P(61)  P(67)  P(71)  P(73)  P(79)  P(83)  P(89)  P(97)
                        P(101) P(103) P(107) P(109) P(113) P(127) P(131) P(137)
                        P(139) P(149) P(151) P(157) P(163) P(167) P(173) P(179)
                        P(181) P(191) P(193) P(197) P(199)
                        default: process_prime_dyn(p, sc, sj, W, vld, rm, k); break;
                    }
                }

                for (uint64_t j = 0; j < W; ++j) {
                    if (vld[j] && rm[j] <= 1) {
                        if (++consec >= k + 1) {
                            uint64_t n = bL + j;
                            if (n > k && n < best.load(std::memory_order_relaxed))
                                if (exact_check(n, k)) {
                                    uint64_t cur = best.load(std::memory_order_relaxed);
                                    while (n < cur && !best.compare_exchange_weak(
                                        cur, n, std::memory_order_relaxed)) {}
                                    found.store(true, std::memory_order_relaxed);
                                }
                        }
                    } else {
                        consec = 0;
                    }
                }
            }
        }
    };

    while (true) {
        found.store(false);
        std::atomic<int> next_chunk{0};
        std::vector<std::thread> threads;
        threads.reserve(g_threads);
        for (unsigned i = 0; i < g_threads; ++i)
            threads.emplace_back(worker, std::ref(next_chunk));
        for (auto& t : threads) t.join();

        if (found.load()) return best.load();
        L_batch += cpb * CHUNK_SIZE;
    }
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------
static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " [-k K] [--start N] [--kmax KMAX] [--threads T]\n";
    std::exit(1);
}

int main(int argc, char** argv) {
    uint64_t k_start = 1, start_L = 1, k_max = 20;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-k" || a == "--k") && i + 1 < argc)
            k_start = std::stoull(argv[++i]);
        else if (a == "--start" && i + 1 < argc)
            start_L = std::stoull(argv[++i]);
        else if (a == "--kmax" && i + 1 < argc)
            k_max = std::stoull(argv[++i]);
        else if (a == "--threads" && i + 1 < argc)
            g_threads = static_cast<unsigned>(std::stoul(argv[++i]));
        else
            usage(argv[0]);
    }

    std::cout << "erdos396-speed | threads=" << g_threads
              << " k=" << k_start << ".." << k_max
              << " start=" << start_L << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    sieve_primes(200'000'000);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Primes: " << g_primes.size() << " in "
              << std::chrono::duration<double>(t1 - t0).count() << "s\n\n";

    for (uint64_t k = k_start; k <= k_max; ++k) {
        auto ts = std::chrono::high_resolution_clock::now();
        uint64_t ans = solve(k, start_L);
        auto te = std::chrono::high_resolution_clock::now();
        double sec = std::chrono::duration<double>(te - ts).count();
        double speed = (ans >= start_L ? ans - start_L + 1 : 1) / sec;

        std::cout << "k=" << std::setw(2) << k
                  << "  n=" << std::setw(15) << ans
                  << "  " << std::fixed << std::setprecision(4) << sec << "s"
                  << "  " << std::setprecision(2) << speed / 1e6 << " M/s\n";
        start_L = ans;
    }
    return 0;
}
