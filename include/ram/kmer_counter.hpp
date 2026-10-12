// Copyright (c) 2024
//
// In-process k-mer counter for the internal-counting path (replaces the
// external FastK count + on-disk table + table load).
//
// Counts canonical k-mers for k in [1, 63]. A k-mer is its packed 2-bit
// canonical encoding (min of the forward/reverse encodings) as a (hi, lo) pair
// of 64-bit words -- exactly the lexicographic min the minimizer engine already
// computes, so the engine can look a k-mer up with the (hi, lo) pair it passes
// today to fastk_get_count(). k <= 32 leaves hi zero and stores only lo.
//
// Singletons are dropped: k-mers occurring fewer than `min_count` times are
// treated as absent (Lookup returns 0). This matches an error-k-mer-dropping
// FastK table and makes whole-genome memory tractable (the error-k-mer
// explosion is overwhelmingly count-1).
//
// Storage: 2^shard_pow2 shards, each an open-addressing (linear probing) table
// of keys and 32-bit counts, 12 bytes per slot (20 for k > 32), counted under
// per-shard locks. With the bloom pre-filter singletons are never stored; without it they
// are stored and dropped in Finalize(), which also shrinks every shard to 3/4
// full. A node-based std::unordered_map took ~75 bytes per k-mer, too much for
// the billions of k-mers of a human genome.

#ifndef RAM_KMER_COUNTER_HPP_
#define RAM_KMER_COUNTER_HPP_

#include <cstddef>
#include <cstdint>
#include <future>
#include <memory>
#include <mutex>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class KmerCounter {
 public:
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

  // Streaming API for overlapping counting with read parsing: CountAsync
  // dispatches one read's counting to the pool and returns immediately (the
  // NucleicAcid must stay alive until Finalize()); Finalize waits for all
  // dispatched work and prunes. Call Finalize() once before any Lookup().
  void CountAsync(const biosoup::NucleicAcid* read);
  void Finalize();

  // O(1) exact lookup of a canonical k-mer given the engine's (hi, lo) words.
  // Returns 0 if the k-mer was pruned or never seen. Lock-free: valid only
  // after Count() has completed (no concurrent writers).
  inline std::uint32_t Lookup(std::uint64_t hi, std::uint64_t lo) const {
    if (!wide_ && hi != 0) {
      return 0;  // k <= 32: every counted k-mer fits in lo
    }
    const std::uint64_t h = KeyHash(hi, lo);
    const Shard& s = shards_[ShardOf(h)];
    if (s.lo.empty()) {
      return 0;
    }
    const std::size_t i = Probe(s, hi, lo, h);
    return IsFree(s, i) ? 0u : s.counts[i];
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

  // splitmix64 finalizer over both key words: every key bit reaches every
  // hash bit. The top shard_pow2_ bits pick the shard, the bits below them
  // the slot.
  static inline std::uint64_t KeyHash(std::uint64_t hi, std::uint64_t lo) {
    std::uint64_t x = lo ^ (hi * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
  }

  // Marks a free slot in `lo` (k <= 32) or `hi` (k > 32). Never part of a
  // canonical k-mer: below k = 32 lo has fewer than 64 bits, for k = 32 all-T
  // is not canonical (all-A is smaller), and for k in [33, 63] hi has at most
  // 62 bits.
  static constexpr std::uint64_t kEmpty = ~0ULL;

  // One shard's table; lo, counts (and hi when wide_) have one entry per slot
  // (none until the first insert). Never full, so a probe always reaches a
  // free slot.
  struct Shard {
    std::mutex m;
    std::vector<std::uint64_t> lo;
    std::vector<std::uint64_t> hi;
    std::vector<std::uint32_t> counts;
    std::size_t size = 0;
  };

  inline std::size_t ShardOf(std::uint64_t h) const {
    return shard_pow2_ == 0 ? 0 : static_cast<std::size_t>(h >> (64 - shard_pow2_));
  }

  inline bool IsFree(const Shard& s, std::size_t i) const {
    return (wide_ ? s.hi[i] : s.lo[i]) == kEmpty;
  }

  // The slot holding (hi, lo) in s, or the free slot where it would go
  // (s must have slots; h = KeyHash(hi, lo)).
  inline std::size_t Probe(const Shard& s, std::uint64_t hi, std::uint64_t lo,
                           std::uint64_t h) const {
    const std::size_t capacity = s.lo.size();
    std::size_t i = static_cast<std::size_t>(
        (static_cast<__uint128_t>(h << shard_pow2_) * capacity) >> 64);
    while (!IsFree(s, i) && !(s.lo[i] == lo && (!wide_ || s.hi[i] == hi))) {
      if (++i == capacity) {
        i = 0;
      }
    }
    return i;
  }

  // Under s.m: the count of (hi, lo), inserted as 0 if absent (grows s first
  // if it would pass 3/4 full).
  std::uint32_t& Slot(Shard& s, std::uint64_t hi, std::uint64_t lo,
                      std::uint64_t h);

  // Rebuild s with `capacity` slots, keeping entries with count >= min_count.
  void Rehash(Shard& s, std::size_t capacity, std::uint32_t min_count);

  // Count one canonical k-mer (from CountRead).
  void Add(std::uint64_t hi, std::uint64_t lo);

  void CountRead(const biosoup::NucleicAcid& read);

  // Drop entries below min_count and shrink every shard to 3/4 full
  // (parallel over shards).
  void Compact();

  // Set the key's bloom bits (lock-free); return true iff all bits were already
  // set (i.e. the key was very likely seen before). Only meaningful if use_bloom_.
  bool BloomTestAndSet(std::uint64_t hi, std::uint64_t lo);

  std::uint32_t k_;
  std::uint32_t min_count_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
  std::uint32_t shard_pow2_;
  bool wide_;  // k > 32: keys need both words
  std::vector<Shard> shards_;
  std::vector<std::future<void>> pending_;  // in-flight CountAsync tasks

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
