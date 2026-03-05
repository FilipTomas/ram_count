#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque types (C++ must not include libfastk.h)
typedef struct Kmer_Table Kmer_Table;
typedef struct Histogram Histogram;

// Loaders
Kmer_Table* ram_fastk_load_kmer_table(const char* name, int cut_off);
Histogram*  ram_fastk_load_histogram(const char* name);

// Query
int64_t ram_fastk_find_kmer(Kmer_Table* table, const char* kseq);
int     ram_fastk_fetch_count(Kmer_Table* table, int64_t idx);

#ifdef __cplusplus
}
#endif