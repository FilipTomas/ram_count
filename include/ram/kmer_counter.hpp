// Copyright (c) 2024
//
// In-process k-mer counter for the internal-counting path (replaces the
// external FastK count + on-disk table + table load).
//
// Counts canonical k-mers for k in [1, 63]. The key is a 128-bit value holding
// the packed 2-bit canonical k-mer (min of the forward/reverse encodings) --
// which is exactly the lexicographic (hi, lo) min the minimizer engine already
// computes, so the engine can look a k-mer up with the (hi, lo) pair it passes
// today to fastk_get_count(). k <= 32 simply leaves the high 64 bits zero.
//
// Singletons are dropped: k-mers occurring fewer than `min_count` times are
// treated as absent (Lookup returns 0). This matches an error-k-mer-dropping
// FastK table and makes whole-genome memory tractable (the error-k-mer
// explosion is overwhelmingly count-1).
//
// Step 1 backend (this file): a sharded, exact hash (std::unordered_map per
// shard) counted under per-shard locks, then pruned below min_count. Exact and
// simple; it holds singletons transiently, so it is meant for chromosome-scale
// validation. A bloom-filtered, bit-packed backend (never storing singletons)
// is intended to slot in behind this same Count()/Lookup() interface for
// whole-genome runs.

#ifndef RAM_KMER_COUNTER_HPP_
#define RAM_KMER_COUNTER_HPP_

#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class KmerCounter {
 public:
  using Key = __uint128_t;

  // k must be in [1, 63]. k-mers seen fewer than `min_count` times are dropped
  // (Lookup returns 0); default 2 ignores singletons.
  //
  // bloom_capacity_hint: if > 0 and min_count >= 2, a bloom pre-filter is used
  // so singletons are never stored in the hash (only their bloom bits). Pass a
  // safe upper bound on the number of distinct k-mers -- total k-mer instances
  // works (distinct <= instances). If 0, all k-mers are stored and pruned at the
  // end (exact, but holds singletons transiently -- not whole-genome viable).
  KmerCounter(
      std::uint32_t k,
      std::shared_ptr<thread_pool::ThreadPool> thread_pool,
      std::uint32_t min_count = 2,
      std::uint64_t bloom_capacity_hint = 0,
      std::uint32_t shard_pow2 = 12);  // 2^12 = 4096 shards

  KmerCounter(const KmerCounter&) = delete;
  KmerCounter& operator=(const KmerCounter&) = delete;

  // Count all canonical k-mers over the reads in [first, last) (multithreaded),
  // then prune entries below min_count. Call once before any Lookup().
  void Count(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last);

  // O(1) exact lookup of a canonical k-mer given the engine's (hi, lo) words.
  // Returns 0 if the k-mer was pruned or never seen. Lock-free: valid only
  // after Count() has completed (no concurrent writers).
  inline std::uint32_t Lookup(std::uint64_t hi, std::uint64_t lo) const {
    const Key key = (static_cast<Key>(hi) << 64) | lo;
    const auto& m = shards_[ShardOf(key)].counts;
    auto it = m.find(key);
    return it == m.end() ? 0u : it->second;
  }

  std::uint32_t k() const { return k_; }
  std::uint64_t num_distinct() const;  // # of retained (>= min_count) k-mers

  // Coverage peak: the count value in [low, high] shared by the most distinct
  // k-mers (low skips the low-count error region). 0 if the table is empty.
  std::uint32_t EstimatePeak(std::uint32_t low = 11,
                             std::uint32_t high = 1000) const;

  // hist[c] = number of distinct k-mers with count c (c in [0, max], capped).
  std::vector<std::uint64_t> Histogram(std::uint32_t max) const;

 private:
  static inline std::uint64_t Mix(std::uint64_t x) {
    x *= 0x9e3779b97f4a7c15ULL;
    x ^= x >> 29;
    return x;
  }
  struct KeyHash {
    std::size_t operator()(Key k) const {
      return static_cast<std::size_t>(
          Mix(static_cast<std::uint64_t>(k)) ^
          Mix(static_cast<std::uint64_t>(k >> 64) + 0x9e3779b97f4a7c15ULL));
    }
  };

  struct Shard {
    std::mutex m;
    std::unordered_map<Key, std::uint32_t, KeyHash> counts;
  };

  inline std::size_t ShardOf(Key key) const {
    return static_cast<std::size_t>(KeyHash{}(key) >> 17) & shard_mask_;
  }

  void CountRead(const std::unique_ptr<biosoup::NucleicAcid>& read);

  // Set the key's bloom bits (lock-free); return true iff all bits were already
  // set (i.e. the key was very likely seen before). Only meaningful if use_bloom_.
  bool BloomTestAndSet(Key key);

  std::uint32_t k_;
  std::uint32_t min_count_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::uint64_t shard_mask_;
  Key mask_;  // (1 << 2k) - 1, masks the packed k-mer to 2k bits
  std::vector<Shard> shards_;

  // Bloom pre-filter (drops singletons without storing them). bloom_ is plain
  // storage accessed atomically via __atomic_* builtins; num bits is a power of
  // two so a mask replaces modulo.
  bool use_bloom_ = false;
  std::vector<std::uint64_t> bloom_;
  std::uint64_t bloom_bit_mask_ = 0;
  std::uint32_t bloom_hashes_ = 4;
};

}  // namespace ram

#endif  // RAM_KMER_COUNTER_HPP_
