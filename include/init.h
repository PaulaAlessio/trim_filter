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
   int trimQ; // NO(0), ALL(1), ENDS(2)
   int trimN; // NO(0), ALL(1), ENDS(2), STRIP(3)
   bool is_fa, is_idx, tree;
   int  minQ, L, Lmer_len;
} Param; 


void printHelpDialog();
Param  get_arg( int argc, char **argv);


#endif
