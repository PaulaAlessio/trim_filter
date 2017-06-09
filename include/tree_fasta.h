//tree_fasta.h 
/* Code that reads in a *fa file and
 * creates a tree containing all the 
 * entries and their reverse complements.*/
#ifndef TREE_FASTA_H
#define TREE_FASTA_H

#include "tree.h"

typedef struct _param_fa{
   long int nlines; 
   long int *entrylen; 
   int linelen;
   int nentries;
   
} Param_fa;


typedef struct _fasta{
   long int  N ;
   unsigned char *seq; 
} Fasta;



int read_fasta(char *filename, Fasta **ptr_fasta);
void free_fasta( Fasta ** ptr_fasta, int n);
void fasta_tree(char * filename, Node* tree, int L);


#endif
