#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SA.h"
#include "Lmer.h"
#include "defines.h"

#define uint64_t unsigned long int 

int cmp2(const void *a, const void *b){
        char* x = &str_input[*(int*)a]; 
        char* y = &str_input[*(int*)b];
        //printf("BEEP %.5s\n", x);
        /* not x>y because I want "aa" before "aab" */
        return strcmp(x,y);
}



int constructSA(int len, int *pos){
      int i ; 
      for (i=0; i<len; i++) pos[i]=i;
      if (len > ((1lu<<32)-1)) 
          return -1; 
      qsort(pos,len,sizeof(int),cmp2);
      return 0;
}


int preindexing(int Lgenome, int Lprefix, int *preindex, int *SA){
  int L_index = 1 << (2*Lprefix);
  int i,j; 
  char prefix[Lprefix];
  int i0=0;
  int cmp; 
  for (i = 0; i < L_index; i++){
     for (j = 0; j < Lprefix; j++){ 
        prefix[Lprefix - j -1] = ( (i>> (2*j))  & 3) + 1;
     }
     int i1 = i0, i2 = Lgenome - 1, i3 = 0;
     while (i2>i1+1){
        i3 =  (i1 + i2) >> 1;
        cmp = strncmp(prefix,str_input + SA[i3],Lprefix);
        if(cmp > 0){
            i1=i3; 
        } else {
            i2=i3;
        }
     }
     if( !strncmp(prefix,str_input + SA[i2],Lprefix) ){
         preindex[i] = i1;
         i0 = i1; 
     }  
     else{
          preindex[i] = -1;
     }
  }
  preindex[L_index] = Lgenome;
  return 1; 

}

int bin_search_read(Fq_read* seq, int Lprefix,int *preindex, int *SA){
      char my_read[seq->L];
      memcpy(my_read,seq->line2,seq->L+1);
      Lmer_sLmer(my_read, seq->L);
      int cmp;
      int flag = 0; 
      while (flag<2){     
         int pos = 0;
         for (int i = 0 ; i < Lprefix; i++){
              if (my_read[i]>'\004') return 0;
              pos += (my_read[i] - 1) << 2*(Lprefix-i -1);
         } 

         //printf("\n%d\n",pos);
         if (preindex[pos] != -1)  {  
            int cc = 1;
            while(preindex[pos+cc]==-1) {
               cc++;
            }
            int i1 = preindex[pos],  i2 = preindex[pos+cc], i3 =0;
            while (i2>i1+1){
              i3 =  (i1 + i2) >> 1;
              cmp = strncmp(&my_read[0],str_input + SA[i3],seq->L);
              if(cmp > 0){
                  i1=i3; 
              } else if (cmp == 0) {
                 return 1;
               } else{
                  i2=i3;
              }
            }
            if (!strncmp(my_read,str_input+ SA[i2],seq->L)){
                return 1;
            }
         } // end if (preindex pos)
         flag++;
         rev_comp2(my_read,seq->L);
      } // end while (check for the reverse complement)   
      return 0;
}

