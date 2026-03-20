#include <stdint.h>
#include <stdlib.h>
#include "libfastk.h"

// Loaders
Kmer_Table* ram_fastk_load_kmer_table(const char* name, int cut_off) {
  return Load_Kmer_Table((char*) name, cut_off);
}

Histogram* ram_fastk_load_histogram(const char* name) {
  return Load_Histogram((char*) name);
}

// K-mer query
int64_t ram_fastk_find_kmer(Kmer_Table* table, const char* kseq) {
  return (int64_t) Find_Kmer(table, (char*) kseq);
}

int ram_fastk_fetch_count(Kmer_Table* table, int64_t idx) {
  return Fetch_Count(table, (int64) idx);
}

// Histogram access
int ram_fastk_hist_low(Histogram* hist) {
  if (hist == NULL) return -1;
  return hist->low;
}

int ram_fastk_hist_high(Histogram* hist) {
  if (hist == NULL) return -1;
  return hist->high;
}

int64_t ram_fastk_hist_value(Histogram* hist, int multiplicity) {
  if (hist == NULL) return -1;
  if (multiplicity < hist->low || multiplicity > hist->high) return 0;
  return (int64_t) hist->hist[multiplicity - hist->low];
}

// Optional cleanup
void ram_fastk_free_kmer_table(Kmer_Table* table) {
  if (table != NULL) {
    Free_Kmer_Table(table);
  }
}

void ram_fastk_free_histogram(Histogram* hist) {
  if (hist != NULL) {
    Free_Histogram(hist);
  }
}