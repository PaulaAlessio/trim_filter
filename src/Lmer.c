#include "Lmer.h"
#define NB 5

#include<stdio.h>
char LT[256]; // global variable. Lookup table
//look up table
void init_map(){
   int i; 
   for(i = 0; i < 256; i++)
      LT[i] = 4; 
   LT['a'] = 0;
   LT['c'] = 1;
   LT['g'] = 2;
   LT['t'] = 3;
   LT['A'] = 0;
   LT['C'] = 1;
   LT['G'] = 2;
   LT['T'] = 3;
}


void Lmer_sLmer(char* Lmer, int L){
   int i;
   for (i = 0; i < L; i++){
      Lmer[i] = LT[(int)Lmer[i]];
   }
}

void rev_comp(char *sLmer, int L){
   char RC[5]= {3,2,1,0,4};
   int c, i, j;
   for (i = 0, j = L-1; i < j; i++, j--) {
      c = RC[(int) sLmer[i]];
      sLmer[i] = RC[(int) sLmer[j]];
      sLmer[j] = c;
   }
}

