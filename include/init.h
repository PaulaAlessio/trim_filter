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
   bool trimQ;
   int trimN; // 0 NO trim, 1 trim all, 2 trim ends
   bool is_fa, is_idx, tree;
   int  minQ, L, Lmer_len;
} Param; 


void printHelpDialog();
Param  get_arg( int argc, char **argv);


#endif
