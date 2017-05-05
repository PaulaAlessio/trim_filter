#ifndef SA_H
#define SA_H
#include "tree_fasta.h" 
#include "fq_read.h"
char *str_input;
int constructSA(int len, int *pos);
int preindexing(int Lgenome, int Lprefix, int *preindex, int *SA);
int bin_search_read(Fq_read *seq, int Lprefix,int *preindex, int *SA);
#endif
