// bloom_maker.h 
// These functions create
#ifndef bloom_maker_H
#define bloom_maker_H

#include "city.h"
#define BITSPERCHAR 8
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
                                             0x40, 0x80 };

typedef struct _bfilter{
   int kmersize, hashNum, kmersizeBytes;
   uint64_t bfsizeBits, bfsizeBytes;
   unsigned char *filter;
} Bfilter;

typedef struct _procs_kmer{
   int kmersize, kmersizeBytes, halfsizeBytes, hangingBases, 
       hasOverhead, hashNum;
   unsigned char *compact;
   uint64_t *hashValues;  
} Procs_kmer;

void init_LUTs();

void init_Bfilter(int kmersize, uint64_t  bfsizeBits, int hashNum, Bfilter *bf);

void free_Bfilter(Bfilter * bf);

void init_procs(Procs_kmer *procs, int kmersize, int hashNum);

void free_procs(Procs_kmer *procs);

int create_bloom_filter(char **fastafiles, int Nfiles , Bfilter *bf);

void save_bloomfilter(char *outputfile, Bfilter *bf);

double is_read_in_filter(unsigned char *read, int L, Procs_kmer *ptr_procs, Bfilter *bf);

#endif
