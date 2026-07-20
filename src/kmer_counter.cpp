// Copyright (c) 2024

#include "ram/kmer_counter.hpp"

#include <future>

namespace ram {

KmerCounter::KmerCounter(
    std::uint32_t k,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::uint32_t min_count,
    std::uint64_t bloom_capacity_hint,
    std::uint32_t shard_pow2)
    : k_(k),
      min_count_(min_count),
      thread_pool_(std::move(thread_pool)),
      shard_mask_((1ULL << shard_pow2) - 1),
      mask_(k >= 64 ? ~static_cast<Key>(0)
                    : ((static_cast<Key>(1) << (2 * k)) - 1)),
      shards_(static_cast<std::size_t>(1) << shard_pow2) {
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

bool KmerCounter::BloomTestAndSet(Key key) {
  // Two 64-bit hashes from the 128-bit key (double hashing); h2 forced odd so
  // (h1 + i*h2) mod 2^b visits distinct positions.
  const std::uint64_t h1 = KeyHash{}(key);
  const std::uint64_t h2 =
      (Mix(static_cast<std::uint64_t>(key >> 64) ^
           (static_cast<std::uint64_t>(key) + 0x100000001b3ULL)) | 1ULL);
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

  Key fwd = 0;
  Key rev = 0;
  const Key mask = mask_;
  const std::uint32_t k = k_;
  const std::uint32_t shift = 2 * (k - 1);

  for (std::uint32_t i = 0; i < read.inflated_len; ++i) {
    const Key c = static_cast<Key>(read.Code(i));          // 0..3
    fwd = ((fwd << 2) | c) & mask;                         // append at the low end
    const Key rc = c ^ 3;                                  // complement
    rev = ((rev >> 2) | (rc << shift)) & mask;             // append complement at the high end

    if (i + 1 >= k) {
      const Key canon = fwd < rev ? fwd : rev;             // canonical = min(fwd, rev)
      Shard& s = shards_[ShardOf(canon)];

      if (use_bloom_) {
        // Fast path: already promoted (seen >= 2 times) -> just increment.
        {
          std::lock_guard<std::mutex> lock(s.m);
          auto it = s.counts.find(canon);
          if (it != s.counts.end()) {
            ++it->second;
            continue;
          }
        }
        // Not promoted yet: promote on the second (bloom-confirmed) sighting;
        // the first sighting only sets bloom bits and is never stored.
        if (BloomTestAndSet(canon)) {
          std::lock_guard<std::mutex> lock(s.m);
          auto res = s.counts.emplace(canon, 2u);
          if (!res.second) {
            ++res.first->second;  // promoted concurrently between find and here
          }
        }
      } else {
        std::lock_guard<std::mutex> lock(s.m);
        ++s.counts[canon];
      }
    }
  }
}

void KmerCounter::CountAsync(const biosoup::NucleicAcid* read) {
  if (thread_pool_ == nullptr) {
    CountRead(*read);
    return;
  }
  pending_.emplace_back(thread_pool_->Submit(
      [this](const biosoup::NucleicAcid* r) -> void { CountRead(*r); }, read));
}

void KmerCounter::Prune() {
  // Drop k-mers below min_count (singletons by default). Safe to run in
  // parallel across shards since each shard is independent.
  if (min_count_ <= 1) {
    return;
  }
  auto prune_shard = [this](Shard& s) {
    for (auto it = s.counts.begin(); it != s.counts.end();) {
      if (it->second < min_count_) {
        it = s.counts.erase(it);
      } else {
        ++it;
      }
    }
  };
  if (thread_pool_ == nullptr) {
    for (auto& s : shards_) prune_shard(s);
  } else {
    std::vector<std::future<void>> futures;
    for (std::size_t i = 0; i < shards_.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&](std::size_t idx) -> void { prune_shard(shards_[idx]); }, i));
    }
    for (auto& f : futures) f.wait();
  }
}

void KmerCounter::Finalize() {
  for (auto& f : pending_) f.wait();
  pending_.clear();
  Prune();
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
    total += s.counts.size();
  }
  return total;
}

std::vector<std::uint64_t> KmerCounter::Histogram(std::uint32_t max) const {
  std::vector<std::uint64_t> hist(static_cast<std::size_t>(max) + 1, 0);
  for (const auto& s : shards_) {
    for (const auto& kv : s.counts) {
      std::uint32_t c = kv.second;
      if (c > max) c = max;
      ++hist[c];
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
