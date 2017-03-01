#include "tree_fasta.h"
#include "defines.h"
#include "fopen_gen.h"
#include <stdlib.h>
#include <string.h>



static Param_fa par;
extern long int alloc_mem; 

static void alloc_param(){
   if (par.nentries % FASTA_ENTRIES == 0){
        par.entrylen = realloc( par.entrylen, 
              sizeof(long int)*(FASTA_ENTRIES  +  par.nentries));
   }  
}

static void init_param(){
   par.nlines = 0;
   par.linelen = 0;
   par.nentries = 0;
   par.entrylen = realloc(par.entrylen,0);
}

static void  get_param_fasta(char *filename){
   FILE *f;
   int offset = 0;  
   char *buffer = (char *) malloc(sizeof(char)*(B_LEN)); 
   int nl_pos = 0; 
   int nc = 0; 
   int newlen = 0;  
   init_param();
   f = fopen_gen(filename,"r");
   if (f == NULL){
      fprintf(stderr, "File %s not found. Exiting program.\n", filename);
      exit(EXIT_FAILURE);
   }
   while ( (newlen = fread( buffer + offset,1, B_LEN - offset,f) ) > 0 ){
      int  j = 0;
      int linelen = 0; 
      newlen += offset;
      //printf("NEWLEN %d\n", newlen);
      while ( j < newlen ){
         switch (buffer[j]){
         // update nlines, linelen, last position after a \n.
         case '\n':
            par.nlines++;
            nl_pos = ++j;
            par.linelen = max(linelen, par.linelen);
            linelen = 0;
            break; 
         //ignore line and annotate number of characters of the entry.
         case '>': 
         case ';':
             alloc_param(); 
             if (par.nlines != 0)
                 par.entrylen[par.nentries++] = nc;
             while(buffer[j] != '\n' &&  j < newlen)
               j++;
             nc = 0  ;  
             break;
         default:
            j++;
            nc++;
            linelen++;
            break;
         }   
      }
      offset = newlen - nl_pos ;
      if(offset > 0){
           memcpy(buffer,buffer + nl_pos, offset);
           nc -= offset;
      }
      nl_pos = 0;
   }
   par.entrylen[par.nentries++] = nc;
   free(buffer);
   fclose(f);
}

// Initialize array of fasta structs (one per entry)
static int  init_fasta(Fasta ** ptr_fasta){
   int i;
   alloc_mem += sizeof(Fasta) * par.nentries; 
   *ptr_fasta = malloc(sizeof(Fasta) * par.nentries);
   for (i = 0; i < par.nentries; i++){
      (*ptr_fasta + i)-> N = par.entrylen[i];
      alloc_mem += sizeof(char) * par.entrylen[i];
      (*ptr_fasta + i)-> seq = malloc(sizeof(char) * par.entrylen[i] );
   }  
   return par.nentries;
}

// ignore header lines
static int ignore_line(char *line){
   int  i = 0;
   while(line[i] != '\n')
      i++;
   return ++i; 
}

int read_fasta(char *filename, Fasta **ptr_fasta){
   
   get_param_fasta(filename);
   int n = init_fasta(ptr_fasta);
   fprintf(stderr, "Reading fasta file: %s.\n- Parameters: \n",filename);
   fprintf(stderr, " * Number of lines: %ld\n", par.nlines );
   fprintf(stderr, " * Number of entries: %d\n", par.nentries);
   fprintf(stderr, " * Length of lines: %d\n",par.linelen);
   fprintf(stderr, " * Length of sequences in entries : [  ");
   int i; 
   for (i = 0; i < n ; i++)
      fprintf(stderr,"%ld ", par.entrylen[i]);
   fprintf(stderr,"].\n");
   // Read the fasta file 
   FILE *f;
   f = fopen_gen(filename,"r");
   if (f == NULL){
      fprintf(stderr, "File %s not found. Exiting program.\n", filename);
      exit(EXIT_FAILURE);
      return -1;
   }
   // Get file size
   fseek(f, 0L, SEEK_END);
   long int  sz = ftell(f);
   fseek(f,0L,SEEK_SET);
   // Allocate memory to read the file in one step
   fprintf(stderr,"- Fasta file size: %ld bytes. \n",sz );
   char*  buffer  =  (char *) malloc(sizeof(char)*sz);
   alloc_mem += sizeof(char)*sz;
   fprintf(stderr,"- Allocating %ld bytes in the buffer. \n",sz );
   mem_usageMB();
   if(sz > 10e8)
      fprintf(stderr,"WARNING: long fasta file. Demanding memory needs.\n");
   long int pos = 0;
   int  linelen =  par.linelen;
   fread(buffer,1,sz,f );
   while(pos < sz){
      for (i = 0; i < n; i++){
         pos += ignore_line(buffer + pos);
         int nfull =  ( (*ptr_fasta + i) -> N )/linelen;
         int remlen =  ( (*ptr_fasta + i) -> N ) % linelen;
         int  j;
         for(j = 0 ; j < nfull; j++){
            memcpy( ( (*ptr_fasta+i)-> seq) + j*linelen, buffer + pos,linelen);
            pos += linelen + 1; 
         } // end read full lines
         // read the remainder
         memcpy( ( (*ptr_fasta+i)-> seq) + j*linelen, buffer + pos,remlen);
         pos += remlen + 1;
      } // end for on the entries
   } // end while on pos
   free(buffer); // free buffer
   alloc_mem -= sizeof(char)*sz;
   fprintf(stderr,"- Deallocating the buffer.\n");
   fprintf(stderr,"- Contents of %s allocated.\n", filename);
   mem_usageMB();
   fclose(f);
   return n; 
}

// free array of fasta structs
void free_fasta( Fasta ** ptr_fasta, int n){
   int i; 
   for (i = 0 ; i < n; i++){
       free((*ptr_fasta + i) -> seq);
       alloc_mem -=sizeof(char)* (*ptr_fasta + i) -> N; 
   }
   free(*ptr_fasta);
   alloc_mem -= sizeof(Fasta) * par.nentries; 

}


// Construct tree from fasta file
void fasta_tree(char *filename,Node *tree, int L){
   Fasta *fasta; 
   int nentries = read_fasta(filename,&fasta);
   int i; 
   fprintf(stderr,"- Constructing  tree ...\n");
   for (i = 0 ; i < nentries; i++){
      construct_tree(tree,(fasta+i)->seq, (fasta+i)-> N, L );
   }
   fprintf(stderr,"- Tree allocated.\n");
   fprintf(stderr,"- Contents of %s deallocated.\n", filename);
   free_fasta(&fasta,nentries);
   fprintf(stderr,"- Fasta contents deallocated.\n");
   mem_usageMB();
}
