#include <stdint.h>
#include "libfastk.h"

// Loaders
Kmer_Table* ram_fastk_load_kmer_table(const char* name, int cut_off) {
  return Load_Kmer_Table((char*)name, cut_off);
}

Histogram* ram_fastk_load_histogram(const char* name) {
  return Load_Histogram((char*)name);
}

// Query
int64_t ram_fastk_find_kmer(Kmer_Table* table, const char* kseq) {
  // FASTK wants char*, but we won't modify the string.
  // If FASTK *does* modify it (unlikely), pass a mutable copy from C++.
  return (int64_t) Find_Kmer(table, (char*)kseq);
}

int ram_fastk_fetch_count(Kmer_Table* table, int64_t idx) {
  // FASTK returns int counts here per your prototype
  return Fetch_Count(table, (int64)idx);
}