#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcat*/
#include "fopen_gen.h"
#include "Lmer.h"
#include "tree.h"
#include "trim_cont.h"
#include "tree_fasta.h" 
#include "defines.h"
#include "fq_read.h"
#include "buf_out.h"
#include "init.h"
#define BUF_FQ 131072
#include <time.h>

extern  char LT[256];
long int alloc_mem;
long int nnodes; 

int main(int argc, char **argv){

   // Get arguments
   Param par = get_arg(argc, argv);

   // Output filenames
   char fq_good[PATH_MAXL]; strcpy(fq_good, par.Oprefix);
   char fq_cont[PATH_MAXL]; strcpy(fq_cont, par.Oprefix);
   char fq_lowq[PATH_MAXL]; strcpy(fq_lowq, par.Oprefix);
   char fq_NNNN[PATH_MAXL]; strcpy(fq_NNNN, par.Oprefix);
   strcat(fq_good,"_good.fq.gz");
   strcat(fq_cont,"_cont.fq.gz");
   strcat(fq_lowq,"_lowq.fq.gz");
   strcat(fq_NNNN,"_NNNN.fq.gz");

   FILE  *fq_in, *f_good, *f_cont, *f_lowq, *f_NNNN;

   int newlen;
   int offset = 0; 
   char buffer[B_LEN + 1];
   int n_cont = 0, n_good = 0, n_lowq = 0, n_trimQ = 0, n_trimN = 0, n_NNNN = 0 ;
   int nreads = 0; 
   int j = 0, k = 0, c1 = 0, c2 = -1;
   bool is_cont;  
   int trimN, trimQ; 
   char char_seq[4*L_LEN]; //string containing one fq read
   int Nchar;   // length of char_seq
   Node* tree =  NULL;

   // time variables
   clock_t start, end;
   double cpu_time_used;
   time_t rawtime;
   struct tm * timeinfo;

   // Start the clock 
   start = clock();
   time ( &rawtime );
   timeinfo = localtime ( &rawtime );
   fprintf(stderr ,"Starting program at: %s", asctime (timeinfo) );

   // Creating tree if flag present
   if(par.tree && !par.is_idx  && par.is_fa){
      fprintf(stderr,"STEP 1: Creating tree.\n");
      init_map();
      tree = new_node_buf();
      fprintf(stderr, "Number of nodes before tree: %ld \n",nnodes);
      fasta_tree(par.Ifa,tree,par.Lmer_len);
      fprintf(stderr, "Number of nodes after tree: %ld \n",nnodes);
   } else if (!par.is_fa){
      fprintf(stderr, "No fasta file given. Not looking for contaminations\n."); 
   }
   else {
      fprintf(stderr,"OPTIONS_WARNING: search method option not implemented yet.\n");
   }

   // Opening files for reading 
   fq_in = fopen_gen(par.Ifq,"r");
   // Open the output files 
   f_good = fopen_gen(fq_good,"w"); 
   f_cont = fopen_gen(fq_cont,"w"); 
   f_lowq = fopen_gen(fq_lowq,"w");
   f_NNNN = fopen_gen(fq_NNNN,"w");
   if( (fq_in == NULL) || (f_good == NULL) ||
       (f_cont == NULL) || (f_lowq == NULL) || (f_NNNN == NULL) ) {
      fprintf(stderr,"Problem opening some of the following files: \n");
      fprintf(stderr," + %s, \n + %s, \n + %s, \n + %s, \n + %s, \nExiting program.\n", 
            par.Ifq,fq_good,fq_cont,fq_lowq, fq_NNNN);
      exit(EXIT_FAILURE);
   }

   // Read the Fq file 
   Fq_read* seq = malloc(sizeof *seq);
   while ( (newlen = fread(buffer+offset,1,B_LEN-offset,fq_in) ) > 0 ){

      newlen += offset+1;
      buffer[newlen] =  '\0';
      for( j = 0 ; buffer[j] != '\0' ; j++){
         c1 = c2 + 1;
         if (buffer[j]== '\n'){
           c2 = j;
           get_sequence(seq,buffer,c1,c2, k);
           if( (k % 4) == 3 ){
              nreads++;
              trimN = trim_sequenceN(seq, par.trimN, par.L>>1);
              if (trimN == 0){
                n_NNNN++;
                Nchar = string_seq(seq, char_seq);
                buffer_output(f_NNNN,char_seq,Nchar,NNNN);
              }else{
                if(trimN == 2) 
                   n_trimN++;
                trimQ = trim_sequenceQ(seq, par.minQ, par.trimQ, par.L>>1 );
                if (trimQ == 0){
                  n_lowq++;
                  Nchar = string_seq(seq, char_seq);
                  buffer_output(f_lowq,char_seq,Nchar,LOWQ);
                } else {
                  if(trimQ == 2)
                     n_trimQ++;
                  is_cont = is_read_in_seq(tree, seq, par.Lmer_len);
                  Nchar = string_seq(seq,char_seq);
                  if(is_cont){
                     n_cont++;
                     buffer_output(f_cont,char_seq,Nchar,CONTAMINATION);
                  } else{
                     n_good++;
                     buffer_output(f_good
                           ,char_seq,Nchar,GOOD);
                  }
                }
              }
              if (nreads % 1000000 == 0) 
                 fprintf(stderr, "  %10d reads have been read.\n",nreads);
           } // end if (k%4 ==3) 
           k++;
         } // end  if \n
     } // end  buffer loop
     offset = newlen - c1 -1;
     if(offset >  -1)
       memcpy(buffer,buffer+c1,offset);
     c2 = -1;
     c1 = 0;
   } // end while 
   // Printing the rest of the buffer outputs
    buffer_output(f_good,NULL,0,GOOD);
    buffer_output(f_NNNN,NULL,0,NNNN);
    buffer_output(f_lowq,NULL,0,LOWQ);
    buffer_output(f_cont,NULL,0,CONTAMINATION);

   // Closing files
   fprintf(stderr, "- Finished reading file.\n");
   fclose(fq_in);
   fclose(f_good);
   fclose(f_cont);  
   fclose(f_lowq);  
   if(par.tree && par.is_fa){ 
      free_all_nodes(); // free memory
   }
   fprintf(stderr, "FINISHING PROGRAM. GENERAL INFO\n");
   fprintf(stderr, "Fastq file had %d reads, of which:\n",nreads); 
   fprintf(stderr, "- %d were discarded due to the presence of N's (-> *_NNNN.fq.gz)\n",n_NNNN); 
   fprintf(stderr, "- %d were discarded due to low quality (-> *_lowq.fq.gz)\n",n_lowq); 
   fprintf(stderr, "- %d were trimmed due to the presence of N's (lowQ and cont filter)\n",n_trimN); 
   fprintf(stderr, "- %d were trimmed due to low quality (-> contamination filter)\n",n_trimQ); 
   fprintf(stderr, "- %d were classified as contaminations (-> *_cont.fq.gz)\n",n_cont); 
   fprintf(stderr, "- %d were classified as good reads (-> *_good.fq.gz)\n",n_good); 
   // Obtaining elapsed time
   end = clock();
   cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC;
   time ( &rawtime );
   timeinfo = localtime ( &rawtime );
   fprintf(stderr ,"Finishing program at: %s", asctime (timeinfo) );
   fprintf(stderr, "Time elapsed: %f s.\n",cpu_time_used);

   return 1; 
}



