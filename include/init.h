//init.h
#ifndef INIT_H
#define INIT_H


#include "defines.h"
/*
 *  Contains structure that reads in the parameters*/

typedef struct param{
   char *Ifq; 
   char *Ifa;
   char *Iidx;
   char *Oprefix;
   int trimQ; // NO(0), FRAC(1), ENDS(2), ENDSFRAC(3), GLOBAL(4)
   int trimN; // NO(0), ALL(1), ENDS(2), STRIP(3)
   bool is_fa, is_idx, tree;
   int  minQ, L, minL, nlowQ, Lmer_len, globleft, globright, percent;
} Param; 


void printHelpDialog();
Param  get_arg( int argc, char **argv);


#endif
