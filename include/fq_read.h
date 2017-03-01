// fq_read.h 
#ifndef FQ_READ_H
#define FQ_READ_H

#include <stdio.h>
#include "tree.h"

#define L_LEN 200
#define ZEROQ 33

typedef struct fq_read{
   char line1[L_LEN], line2[L_LEN], line3[L_LEN], line4[L_LEN];
   int L;
} Fq_read;

void get_sequence(Fq_read* seq, char* buffer,int c1, int c2, int k);


int string_seq( Fq_read *seq, char *char_seq );

#endif
