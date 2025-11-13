#include "biosoup/nucleic_acid.hpp"
#include <cstdint>
#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>

struct ReadRec {
  std::unique_ptr<biosoup::NucleicAcid> seq;
  std::unordered_map<std::uint64_t, std::uint64_t> kmer_counts;
};