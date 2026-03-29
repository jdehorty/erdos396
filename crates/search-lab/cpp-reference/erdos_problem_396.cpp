#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <thread>
#include <chrono>
#include <iomanip>
#include <bit>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <cstdlib>

unsigned int NUM_THREADS = std::thread::hardware_concurrency();

struct PrimeData
{
    uint64_t magic;
    uint64_t inv_p;
    uint64_t limit;
    uint32_t p;
    uint8_t shift;
};

struct BucketItem
{
    uint32_t p_idx;
    uint32_t offset;
};

struct FastBucket
{
    BucketItem* data;
    uint32_t count;
    uint32_t cap;

    FastBucket()
    {
        cap = 16384;
        data = static_cast<BucketItem*>(std::malloc(cap * sizeof(BucketItem)));
        count = 0;
    }

    ~FastBucket() { std::free(data); }

    FastBucket(const FastBucket&) = delete;
    FastBucket& operator=(const FastBucket&) = delete;

    inline void clear() { count = 0; }

    inline void push_back(uint32_t p_idx, uint32_t offset)
    {
        if (count == cap) [[unlikely]]
        {
            cap *= 2;
            data = static_cast<BucketItem*>(std::realloc(data, cap * sizeof(BucketItem)));
        }
        data[count].p_idx = p_idx;
        data[count].offset = offset;
        ++count;
    }
};

std::vector<PrimeData> primes;

void get_primes(uint32_t limit)
{
    std::vector<bool> sieve(limit + 1, true);
    sieve[0] = false;
    sieve[1] = false;
    for (uint64_t p = 2; p * p <= limit; ++p)
    {
        if (sieve[p])
        {
            for (uint64_t i = p * p; i <= limit; i += p)
                sieve[i] = false;
        }
    }
    for (uint32_t p = 2; p <= limit; ++p)
    {
        if (sieve[p])
        {
            uint8_t shift = std::bit_width(p) - 1;
            unsigned __int128 magic_128 = ((unsigned __int128)1 << (64 + shift)) / p + 1;

            uint64_t inv = 0;
            if (p % 2 != 0)
            {
                inv = p;
                for (int i = 0; i < 5; ++i)
                    inv *= 2 - p * inv;
            }
            primes.push_back({(uint64_t)magic_128, inv, UINT64_MAX / p, p, shift});
        }
    }
}

bool exact_check(uint64_t n, uint64_t k)
{
    uint32_t nu_2_prod = std::popcount(n - k - 1) - std::popcount(n) + k + 1;
    if (std::popcount(n) < nu_2_prod) [[unlikely]]
        return false;

    for (const auto& pd : primes)
    {
        uint32_t p = pd.p;
        if (p == 2)
            continue;
        if (p > 2 * k)
            break;

        uint64_t nu_prod = 0, nu_comb = 0, power = p;
        while (true)
        {
            uint64_t v_n = n / power;
            nu_prod += v_n - (n - k - 1) / power;
            nu_comb += (2 * n) / power - 2 * v_n;
            if (power > (2 * n) / p)
                break;
            power *= p;
        }
        if (nu_prod > nu_comb)
            return false;
    }

    for (uint64_t i = 0; i <= k; ++i)
    {
        uint64_t temp = n - i;
        temp >>= std::countr_zero(temp);
        for (const auto& pd : primes)
        {
            uint32_t p = pd.p;
            if (p == 2)
                continue;

            if (p <= 2 * k)
            {
                uint64_t q = temp * pd.inv_p;
                while (q <= pd.limit)
                {
                    temp = q;
                    q = temp * pd.inv_p;
                }
                continue;
            }
            if (p * p > temp)
            {
                if (temp > 2 * k)
                {
                    uint64_t p_val = temp;
                    uint64_t nu_prod = 1;
                    uint64_t t2 = (n - i) / p_val;
                    while (t2 > 0 && t2 % p_val == 0)
                    {
                        nu_prod++;
                        t2 /= p_val;
                    }

                    uint64_t nu_comb = 0;
                    uint64_t power = p_val;
                    while (true)
                    {
                        uint64_t v_n = n / power;
                        nu_comb += (2 * n) / power - 2 * v_n;
                        if (power > (2 * n) / p_val)
                            break;
                        power *= p_val;
                    }
                    if (nu_prod > nu_comb)
                        return false;
                }
                break;
            }
            uint64_t q = temp * pd.inv_p;
            if (q <= pd.limit)
            {
                uint64_t nu_prod = 0;
                while (q <= pd.limit)
                {
                    nu_prod++;
                    temp = q;
                    q = temp * pd.inv_p;
                }

                uint64_t nu_comb = 0;
                uint64_t power = p;
                while (true)
                {
                    uint64_t v_n = n / power;
                    nu_comb += (2 * n) / power - 2 * v_n;
                    if (power > (2 * n) / p)
                        break;
                    power *= p;
                }
                if (nu_prod > nu_comb)
                    return false;
            }
        }
    }
    return true;
}

template<uint32_t p>
inline void process_prime_p(uint32_t& start_j, uint32_t W_block, uint64_t* rem)
{
    constexpr uint64_t inv_p = []() {
        uint64_t inv = p;
        for (int i = 0; i < 5; ++i)
            inv *= 2 - p * inv;
        return inv;
    }();
    constexpr uint64_t limit = UINT64_MAX / p;

    uint32_t j = start_j;
    for (; j < W_block; j += p)
    {
        uint64_t temp = rem[j] * inv_p;
        uint64_t q = temp * inv_p;
        if (q <= limit) [[unlikely]]
        {
            temp = q;
            while (true)
            {
                q = temp * inv_p;
                if (q > limit)
                    break;
                temp = q;
            }
        }
        rem[j] = temp;
    }
    start_j = j - W_block;
}

inline void process_p_dyn(uint32_t p, uint64_t inv_p, uint64_t limit, uint32_t& start_j, uint32_t W_block, uint64_t* rem)
{
    uint32_t j = start_j;
    for (; j < W_block; j += p)
    {
        uint64_t temp = rem[j] * inv_p;
        uint64_t q = temp * inv_p;
        if (q <= limit) [[unlikely]]
        {
            temp = q;
            while (true)
            {
                q = temp * inv_p;
                if (q > limit)
                    break;
                temp = q;
            }
        }
        rem[j] = temp;
    }
    start_j = j - W_block;
}

struct SolveResult {
    uint64_t min_witness;
    uint64_t witness_count;
};

SolveResult solve(uint64_t k, uint64_t start_L, uint64_t end_L = UINT64_MAX, bool bench_mode = false)
{
    const uint64_t CHUNK_SIZE = 1048576;
    std::atomic<uint64_t> current_chunk{0};
    std::atomic<uint64_t> global_min_n{static_cast<uint64_t>(-1)};
    std::atomic<uint64_t> witness_count{0};

    auto worker = [&]() {
        const uint32_t OPT_BLOCK_SIZE = 32768;
        const uint32_t BLOCK_SHIFT = 15;
        const uint32_t BLOCK_MASK = 32767;

        std::vector<uint64_t> rem(OPT_BLOCK_SIZE + 32);
        std::vector<uint32_t> prime_offsets;

        FastBucket buckets[33];

        while (true)
        {
            uint64_t chunk_id = current_chunk.fetch_add(1, std::memory_order_relaxed);
            uint64_t L_chunk = start_L + chunk_id * CHUNK_SIZE;
            uint64_t R_chunk = L_chunk + CHUNK_SIZE + k - 1;

            if (L_chunk >= end_L)
                break;
            if (!bench_mode && L_chunk > global_min_n.load(std::memory_order_relaxed))
                break;

            uint32_t CHUNK_W = (uint32_t)(R_chunk - L_chunk + 1);
            uint64_t max_p = std::max<uint64_t>(std::sqrt(2 * R_chunk) + 1, 2 * k);
            auto it = std::upper_bound(primes.begin(), primes.end(), max_p,[](uint64_t val, const PrimeData& pd) { return val < pd.p; });
            size_t chunk_total_primes = std::distance(primes.begin(), it);

            if (prime_offsets.size() < chunk_total_primes)
                prime_offsets.resize(chunk_total_primes);

            auto it_large = std::upper_bound(primes.begin(), primes.begin() + chunk_total_primes, OPT_BLOCK_SIZE,[](uint32_t val, const PrimeData& pd) { return val < pd.p; });
            size_t first_large_prime_idx = std::distance(primes.begin(), it_large);

            for (size_t idx = 1; idx < chunk_total_primes; ++idx)
            {
                uint32_t p = primes[idx].p;
                uint64_t magic = primes[idx].magic;
                uint8_t shift = primes[idx].shift;

                uint64_t num = L_chunk + p - 1;
                uint64_t start_c = (uint64_t)((((unsigned __int128)num * magic) >> 64) >> shift);
                prime_offsets[idx] = (uint32_t)(start_c * p - L_chunk);
            }

            for (int i = 0; i < 33; ++i) buckets[i].clear();

            for (size_t idx = first_large_prime_idx; idx < chunk_total_primes; ++idx)
            {
                uint64_t start_j = prime_offsets[idx];
                uint32_t p = primes[idx].p;
                while (start_j < CHUNK_W)
                {
                    buckets[start_j >> BLOCK_SHIFT].push_back((uint32_t)idx, (uint32_t)(start_j & BLOCK_MASK));
                    start_j += p;
                }
            }

            uint32_t overlap = 0;
            uint32_t j = 0;

            for (uint32_t block_start = 0; block_start < CHUNK_W; block_start += OPT_BLOCK_SIZE)
            {
                uint32_t block_idx = block_start >> BLOCK_SHIFT;
                uint64_t block_L = L_chunk + block_start;
                uint64_t block_R = std::min(R_chunk, block_L + OPT_BLOCK_SIZE - 1);
                uint32_t num_new = (uint32_t)(block_R - block_L + 1);

                uint64_t x = block_L;
                for (uint32_t idx_new = 0; idx_new < num_new; ++idx_new, ++x)
                    rem[overlap + idx_new] = x >> std::countr_zero(x);

                uint64_t* rem_ptr = rem.data() + overlap;

                for (size_t idx = 1; idx < first_large_prime_idx; ++idx)
                {
                    uint32_t p = primes[idx].p;
                    uint32_t start_j = prime_offsets[idx];

                    if (start_j < num_new)
                    {
                        uint64_t inv_p = primes[idx].inv_p;
                        uint64_t limit = primes[idx].limit;

                        #define PROCESS_PRIME(p_val) \
                            case p_val: process_prime_p<p_val>(start_j, num_new, rem_ptr); break;

                        switch (p)
                        {
                            PROCESS_PRIME(3) PROCESS_PRIME(5) PROCESS_PRIME(7) PROCESS_PRIME(11)
                            PROCESS_PRIME(13) PROCESS_PRIME(17) PROCESS_PRIME(19) PROCESS_PRIME(23)
                            PROCESS_PRIME(29) PROCESS_PRIME(31) PROCESS_PRIME(37) PROCESS_PRIME(41)
                            PROCESS_PRIME(43) PROCESS_PRIME(47) PROCESS_PRIME(53) PROCESS_PRIME(59)
                            PROCESS_PRIME(61) PROCESS_PRIME(67) PROCESS_PRIME(71) PROCESS_PRIME(73)
                            PROCESS_PRIME(79) PROCESS_PRIME(83) PROCESS_PRIME(89) PROCESS_PRIME(97)
                            default: process_p_dyn(p, inv_p, limit, start_j, num_new, rem_ptr);
                        }
                        #undef PROCESS_PRIME
                        prime_offsets[idx] = start_j;
                    }
                    else
                    {
                        prime_offsets[idx] = start_j - num_new;
                    }
                }

                uint32_t b_count = buckets[block_idx].count;
                BucketItem* b_items = buckets[block_idx].data;

                for (uint32_t i = 0; i < b_count; ++i)
                {
                    if (i + 8 < b_count) [[likely]]
                        __builtin_prefetch(&primes[b_items[i + 8].p_idx], 0, 1);

                    uint32_t p_idx = b_items[i].p_idx;
                    uint32_t offset = b_items[i].offset;

                    uint64_t inv_p = primes[p_idx].inv_p;
                    uint64_t limit = primes[p_idx].limit;

                    uint64_t temp = rem_ptr[offset] * inv_p;
                    uint64_t q = temp * inv_p;
                    if (q <= limit) [[unlikely]]
                    {
                        temp = q;
                        while (true)
                        {
                            q = temp * inv_p;
                            if (q > limit)
                                break;
                            temp = q;
                        }
                    }
                    rem_ptr[offset] = temp;
                }

                uint32_t W_search = overlap + num_new;
                while (j + k < W_search)
                {
                    int i = k;
                    while (i >= 0 && rem[j + i] <= 1)
                        --i;

                    if (i < 0)
                    {
                        uint64_t n = block_L - overlap + j + k;
                        if (bench_mode)
                        {
                            if (n > k && exact_check(n, k))
                            {
                                witness_count.fetch_add(1, std::memory_order_relaxed);
                                uint64_t current = global_min_n.load(std::memory_order_relaxed);
                                while (n < current && !global_min_n.compare_exchange_weak(current, n, std::memory_order_relaxed)) {}
                            }
                        }
                        else if (n > k && n < global_min_n.load(std::memory_order_relaxed))
                        {
                            if (exact_check(n, k))
                            {
                                uint64_t current = global_min_n.load(std::memory_order_relaxed);
                                while (n < current && !global_min_n.compare_exchange_weak(current, n, std::memory_order_relaxed)) {}
                            }
                        }
                        ++j;
                    }
                    else
                    {
                        j += i + 1;
                    }
                }

                if (W_search >= k)
                {
                    for (uint32_t i = 0; i < k; ++i)
                        rem[i] = rem[W_search - k + i];
                    j -= (W_search - k);
                    overlap = k;
                }
                else
                {
                    overlap = W_search;
                }
            }

            if (!bench_mode)
            {
                static std::atomic<bool> is_writing{false};
                if (chunk_id > NUM_THREADS && chunk_id % 1000 == 0)
                {
                    bool expected = false;
                    if (is_writing.compare_exchange_strong(expected, true))
                    {
                        uint64_t safe_L = start_L + (chunk_id - NUM_THREADS) * CHUNK_SIZE;
                        std::ofstream fout("checkpoint-396.tmp");
                        if (fout)
                        {
                            fout << k << " " << safe_L << "\n";
                            fout.close();
                            std::error_code ec;
                            std::filesystem::rename("checkpoint-396.tmp", "checkpoint-396.txt", ec);
                        }
                        is_writing.store(false, std::memory_order_release);
                    }
                }
            }
        }
    };

    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < NUM_THREADS; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads)
        t.join();

    return {global_min_n.load(), witness_count.load()};
}

int main(int argc, char** argv)
{
    uint64_t start_k = 0, start_L = 0, end_L = UINT64_MAX, k_max = 20;
    bool bench_mode = false;

    // Parse CLI args: -k K, --start N, --end N, --kmax K, --threads T, --bench SECS
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if ((a == "-k" || a == "--k") && i + 1 < argc)
            start_k = std::stoull(argv[++i]);
        else if (a == "--start" && i + 1 < argc)
            start_L = std::stoull(argv[++i]);
        else if (a == "--end" && i + 1 < argc)
            end_L = std::stoull(argv[++i]);
        else if (a == "--kmax" && i + 1 < argc)
            k_max = std::stoull(argv[++i]);
        else if (a == "--threads" && i + 1 < argc)
            NUM_THREADS = static_cast<unsigned int>(std::stoull(argv[++i]));
        else if (a == "--bench" && i + 1 < argc)
        {
            ++i; // consume SECS (unused — bench mode just ignores witnesses)
            bench_mode = true;
        }
        else
        {
            std::cerr << "Usage: " << argv[0]
                      << " [-k K] [--start N] [--end N] [--kmax K] [--threads T] [--bench SECS]\n";
            return 1;
        }
    }

    // Fall back to checkpoint if no CLI args given
    if (start_k == 0)
    {
        start_k = 1;
        start_L = 1;
        if (std::filesystem::exists("checkpoint-396.txt"))
        {
            std::ifstream fin("checkpoint-396.txt");
            if (!(fin >> start_k >> start_L))
            {
                start_k = 1;
                start_L = 1;
            }
            else
            {
                std::cout << "--> Resuming from checkpoint: k = " << start_k
                          << ", L_batch = " << start_L << "\n\n";
            }
        }
    }

    std::cout << "# erdos396-ref | threads=" << NUM_THREADS
              << " k=" << start_k << ".." << k_max
              << " start=" << start_L
              << " end=" << (end_L == UINT64_MAX ? "inf" : std::to_string(end_L))
              << "\n";

    auto start_primes = std::chrono::high_resolution_clock::now();
    get_primes(200000000);
    auto end_primes = std::chrono::high_resolution_clock::now();
    std::cout << "# primes: " << primes.size() << " in "
              << std::fixed << std::setprecision(3)
              << std::chrono::duration<double>(end_primes - start_primes).count()
              << "s\n";

    for (uint64_t k = start_k; k <= k_max; ++k)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = solve(k, start_L, end_L, bench_mode);
        auto end = std::chrono::high_resolution_clock::now();

        uint64_t ans = result.min_witness;
        uint64_t witnesses = result.witness_count;
        std::chrono::duration<double> elapsed = end - start;
        double seconds = elapsed.count();

        uint64_t candidates_checked;
        std::string witness_str;
        if (bench_mode && end_L < UINT64_MAX)
        {
            candidates_checked = end_L - start_L;
            witness_str = "bench";
        }
        else if (ans < UINT64_MAX && ans >= start_L)
        {
            candidates_checked = ans - start_L + 1;
            witness_str = std::to_string(ans);
        }
        else if (end_L < UINT64_MAX)
        {
            candidates_checked = end_L - start_L;
            witness_str = "none";
        }
        else
        {
            candidates_checked = 1;
            witness_str = "none";
        }
        double speed = candidates_checked / seconds;

        std::ostringstream oss;
        oss << "R\t" << k << "\t" << witness_str
            << "\t" << std::fixed << std::setprecision(4) << seconds
            << "\t" << candidates_checked
            << "\t" << std::fixed << std::setprecision(2) << (speed / 1e6)
            << "\t" << witnesses << "\n";

        std::string output_str = oss.str();
        std::cout << output_str;

        if (!bench_mode)
        {
            std::ofstream results_file("results-396.txt", std::ios::app);
            if (results_file)
                results_file << output_str;

            if (ans < UINT64_MAX)
                start_L = ans;
            else
                break;

            std::ofstream fout("checkpoint-396.tmp");
            if (fout)
            {
                fout << (k + 1) << " " << start_L << "\n";
                fout.close();
                std::error_code ec;
                std::filesystem::rename("checkpoint-396.tmp", "checkpoint-396.txt", ec);
            }
        }
    }

    std::error_code ec;
    std::filesystem::remove("checkpoint-396.txt", ec);

    return 0;
}