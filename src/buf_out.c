#include <stdio.h>
#include <string.h>
#include "buf_out.h"
#include "defines.h"



void buffer_output(FILE *fout, const char*a, const int len, const int fd_i){
   //defined static so that it doesn't get inizialized every time
   static char buf[4][B_LEN]; 
   static int count[4] = {0,0,0,0} ;

   // empties the buffer if there is no "a" or is buffer is full or if len == 0
   if (count[fd_i] + len >= B_LEN || a==NULL|| strlen(a)!=len || len==0){
      fwrite(buf[fd_i], 1, count[fd_i], fout);
      count[fd_i] = 0 ; 
   }
   memcpy(buf[fd_i]+count[fd_i],a, len);
   count[fd_i]+= len; 
}

void write_summary(int nreads, int n_good, int n_NNNN, int n_lowq, 
      int n_trimN, int n_trimQ, int n_cont, char *summary){

      FILE *f; 
      f = fopen(summary,"wb");
      fwrite(&nreads,sizeof(int),1,f);
      fwrite(&n_good,sizeof(int),1,f);
      fwrite(&n_NNNN,sizeof(int),1,f);
      fwrite(&n_lowq,sizeof(int),1,f);
      fwrite(&n_trimN,sizeof(int),1,f);
      fwrite(&n_trimQ,sizeof(int),1,f);
      fwrite(&n_cont,sizeof(int),1,f);
      fclose(f); 

} 
