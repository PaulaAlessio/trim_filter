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


