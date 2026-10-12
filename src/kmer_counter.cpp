// Copyright (c) 2024

#include "ram/kmer_counter.hpp"

#include <algorithm>
#include <future>
#include <stdexcept>
#include <string>

namespace ram {

constexpr std::uint64_t KmerCounter::kEmpty;

KmerCounter::KmerCounter(
    std::uint32_t k,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::uint32_t min_count,
    std::uint64_t bloom_capacity_hint,
    std::uint32_t shard_pow2)
    : k_(k),
      min_count_(min_count),
      thread_pool_(std::move(thread_pool)),
      shard_pow2_(shard_pow2),
      wide_(k > 32),
      shards_(static_cast<std::size_t>(1) << shard_pow2) {
  if (k < 1 || k > 63) {
    throw std::invalid_argument(
        "[ram::KmerCounter] error: k must be in [1, 63], got " + std::to_string(k));
  }
  if (shard_pow2 > 32) {
    throw std::invalid_argument(
        "[ram::KmerCounter] error: shard_pow2 must be at most 32");
  }
  if (bloom_capacity_hint > 0 && min_count_ >= 2) {
    // num bits = next power of two >= hint (>= 1 bit per expected distinct
    // k-mer; using total instances as the hint gives many bits per distinct).
    std::uint64_t num_bits = 64;
    while (num_bits < bloom_capacity_hint) {
      num_bits <<= 1;
    }
    bloom_.assign(static_cast<std::size_t>(num_bits / 64), 0);
    bloom_bit_mask_ = num_bits - 1;
    bloom_hashes_ = 4;
    use_bloom_ = true;
  }
}

bool KmerCounter::BloomTestAndSet(std::uint64_t hi, std::uint64_t lo) {
  // Two 64-bit hashes from the key (double hashing); h2 forced odd so
  // (h1 + i*h2) mod 2^b visits distinct positions.
  const std::uint64_t h1 = Mix(lo) ^ Mix(hi + 0x9e3779b97f4a7c15ULL);
  const std::uint64_t h2 = Mix(hi ^ (lo + 0x100000001b3ULL)) | 1ULL;
  bool all_already_set = true;
  for (std::uint32_t i = 0; i < bloom_hashes_; ++i) {
    const std::uint64_t bit = (h1 + static_cast<std::uint64_t>(i) * h2) & bloom_bit_mask_;
    std::uint64_t* word = &bloom_[static_cast<std::size_t>(bit >> 6)];
    const std::uint64_t m = 1ULL << (bit & 63);
    const std::uint64_t old = __atomic_fetch_or(word, m, __ATOMIC_RELAXED);
    if ((old & m) == 0) {
      all_already_set = false;  // this bit was set by us now -> first time
    }
  }
  return all_already_set;
}

void KmerCounter::CountRead(const biosoup::NucleicAcid& read) {
  if (read.inflated_len < k_) {
    return;
  }

  const std::uint32_t k = k_;
  const std::uint32_t shift = 2 * (k - 1);

  if (!wide_) {
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    const std::uint64_t mask = k == 32 ? ~0ULL : (1ULL << (2 * k)) - 1;
    for (std::uint32_t i = 0; i < read.inflated_len; ++i) {
      const std::uint64_t c = read.Code(i);                // 0..3
      fwd = ((fwd << 2) | c) & mask;                       // append at the low end
      rev = ((rev >> 2) | ((c ^ 3) << shift)) & mask;      // complement at the high end
      if (i + 1 >= k) {
        Add(0, fwd < rev ? fwd : rev);                     // canonical = min(fwd, rev)
      }
    }
    return;
  }

  __uint128_t fwd = 0;
  __uint128_t rev = 0;
  const __uint128_t mask = (static_cast<__uint128_t>(1) << (2 * k)) - 1;
  for (std::uint32_t i = 0; i < read.inflated_len; ++i) {
    const __uint128_t c = read.Code(i);
    fwd = ((fwd << 2) | c) & mask;
    rev = ((rev >> 2) | ((c ^ 3) << shift)) & mask;
    if (i + 1 >= k) {
      const __uint128_t canon = fwd < rev ? fwd : rev;
      Add(static_cast<std::uint64_t>(canon >> 64), static_cast<std::uint64_t>(canon));
    }
  }
}

void KmerCounter::Add(std::uint64_t hi, std::uint64_t lo) {
  const std::uint64_t h = KeyHash(hi, lo);
  Shard& s = shards_[ShardOf(h)];

  if (use_bloom_) {
    // Fast path: already promoted (seen >= 2 times) -> just increment.
    {
      std::lock_guard<std::mutex> lock(s.m);
      if (!s.lo.empty()) {
        const std::size_t i = Probe(s, hi, lo, h);
        if (!IsFree(s, i)) {
          ++s.counts[i];
          return;
        }
      }
    }
    // Not promoted yet: promote on the second (bloom-confirmed) sighting;
    // the first sighting only sets bloom bits and is never stored.
    if (BloomTestAndSet(hi, lo)) {
      std::lock_guard<std::mutex> lock(s.m);
      std::uint32_t& count = Slot(s, hi, lo, h);
      // nonzero: promoted concurrently between the two locks
      count = count == 0 ? 2u : count + 1;
    }
  } else {
    std::lock_guard<std::mutex> lock(s.m);
    ++Slot(s, hi, lo, h);
  }
}

std::uint32_t& KmerCounter::Slot(Shard& s, std::uint64_t hi, std::uint64_t lo,
                                 std::uint64_t h) {
  if ((s.size + 1) * 4 > s.lo.size() * 3) {
    Rehash(s, std::max<std::size_t>(1024, 2 * s.lo.size()), 0);
  }
  const std::size_t i = Probe(s, hi, lo, h);
  if (IsFree(s, i)) {
    s.lo[i] = lo;
    if (wide_) {
      s.hi[i] = hi;
    }
    s.counts[i] = 0;
    ++s.size;
  }
  return s.counts[i];
}

void KmerCounter::Rehash(Shard& s, std::size_t capacity,
                         std::uint32_t min_count) {
  Shard next;
  next.lo.assign(capacity, kEmpty);
  if (wide_) {
    next.hi.assign(capacity, kEmpty);
  }
  next.counts.assign(capacity, 0);
  for (std::size_t i = 0; i < s.lo.size(); ++i) {
    if (IsFree(s, i) || s.counts[i] < min_count) {
      continue;
    }
    const std::uint64_t hi = wide_ ? s.hi[i] : 0;
    const std::size_t j = Probe(next, hi, s.lo[i], KeyHash(hi, s.lo[i]));
    next.lo[j] = s.lo[i];
    if (wide_) {
      next.hi[j] = hi;
    }
    next.counts[j] = s.counts[i];
    ++next.size;
  }
  s.lo.swap(next.lo);
  s.hi.swap(next.hi);
  s.counts.swap(next.counts);
  s.size = next.size;
}

void KmerCounter::CountAsync(const biosoup::NucleicAcid* read) {
  if (thread_pool_ == nullptr) {
    CountRead(*read);
    return;
  }
  pending_.emplace_back(thread_pool_->Submit(
      [this](const biosoup::NucleicAcid* r) -> void { CountRead(*r); }, read));
}

void KmerCounter::Compact() {
  // Shards are independent, so this runs in parallel across them.
  auto compact_shard = [this](Shard& s) {
    std::size_t kept = 0;
    for (std::size_t i = 0; i < s.lo.size(); ++i) {
      if (!IsFree(s, i) && s.counts[i] >= min_count_) {
        ++kept;
      }
    }
    if (kept == 0) {
      std::vector<std::uint64_t>().swap(s.lo);
      std::vector<std::uint64_t>().swap(s.hi);
      std::vector<std::uint32_t>().swap(s.counts);
      s.size = 0;
      return;
    }
    Rehash(s, kept * 4 / 3 + 1, min_count_);
  };
  if (thread_pool_ == nullptr) {
    for (auto& s : shards_) compact_shard(s);
  } else {
    std::vector<std::future<void>> futures;
    for (std::size_t i = 0; i < shards_.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&](std::size_t idx) -> void { compact_shard(shards_[idx]); }, i));
    }
    for (auto& f : futures) f.wait();
  }
}

void KmerCounter::Finalize() {
  for (auto& f : pending_) f.wait();
  pending_.clear();
  Compact();
}

void KmerCounter::Count(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last) {
  for (auto it = first; it != last; ++it) {
    CountAsync(it->get());
  }
  Finalize();
}

std::uint64_t KmerCounter::num_distinct() const {
  std::uint64_t total = 0;
  for (const auto& s : shards_) {
    total += s.size;
  }
  return total;
}

std::vector<std::uint64_t> KmerCounter::Histogram(std::uint32_t max) const {
  // Every stride-th shard, from shard `first`; hash-map nodes are scattered,
  // so one thread walking all of them takes minutes on a human chromosome.
  auto shard_histogram = [&](std::size_t first, std::size_t stride) {
    std::vector<std::uint64_t> hist(static_cast<std::size_t>(max) + 1, 0);
    for (std::size_t i = first; i < shards_.size(); i += stride) {
      const Shard& s = shards_[i];
      for (std::size_t j = 0; j < s.lo.size(); ++j) {
        if (!IsFree(s, j)) {
          ++hist[std::min(s.counts[j], max)];
        }
      }
    }
    return hist;
  };
  if (thread_pool_ == nullptr) {
    return shard_histogram(0, 1);
  }
  const std::size_t stride =
      std::min<std::size_t>(thread_pool_->num_threads(), shards_.size());
  std::vector<std::future<std::vector<std::uint64_t>>> futures;
  for (std::size_t i = 0; i < stride; ++i) {
    futures.emplace_back(thread_pool_->Submit(shard_histogram, i, stride));
  }
  std::vector<std::uint64_t> hist(static_cast<std::size_t>(max) + 1, 0);
  for (auto& f : futures) {
    const std::vector<std::uint64_t> part = f.get();
    for (std::size_t c = 0; c < hist.size(); ++c) {
      hist[c] += part[c];
    }
  }
  return hist;
}

std::uint32_t KmerCounter::EstimatePeak(std::uint32_t low,
                                        std::uint32_t high) const {
  const std::vector<std::uint64_t> hist = Histogram(high);
  std::uint32_t best = 0;
  std::uint64_t best_v = 0;
  for (std::uint32_t c = low; c <= high; ++c) {
    if (hist[c] > best_v) {
      best_v = hist[c];
      best = c;
    }
  }
  return best;
}

}  // namespace ram
