//trim_cont.h 
#ifndef TRIM_CONT_H
#define TRIM_CONT_H

#include "fq_read.h"
#include "defines.h"

/* returns 0 if no N's found, 1 if N's found*/
int no_N(Fq_read *seq);

/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int Nfree_Lmer(Fq_read *seq, int minL);

/* returns 1 if no N's found, 2 if trimmed*/
int Ntrim_ends(Fq_read *seq, int minL);

/* 0 if sequence contains low quality nucleotides, 1 otherwise*/
int no_lowQ(Fq_read *seq, int minQ); 

/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int Qtrim_ends(Fq_read *seq, int minQ, int minL);

/* Checking N  base callings*/
/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceN(Fq_read *seq, int trimN, int minL);

/* Checking read quality */
/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceQ(Fq_read *seq, int minQ, int trimQ, int minL);

/* Checks if a read is in the tree: true if found, false otherwise*/
bool is_read_in_seq( Node *tree, Fq_read *seq, int L);

#endif