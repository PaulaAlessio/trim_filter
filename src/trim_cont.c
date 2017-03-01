#include <string.h>
#include "trim_cont.h"
#include "Lmer.h"

/* returns 0 if no N's found, 1 if N's found*/
int no_N(Fq_read *seq){
   char Lmer[seq-> L];
   memcpy(Lmer,seq->line2,seq -> L);
   Lmer_sLmer(Lmer,seq->L);
   for(int i = 0; i < seq->L; i++){
      if ( Lmer[i]  > '\003' ){
         return 0;
      }
   }
   return 1;
}


/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int Nfree_Lmer(Fq_read *seq, int minL){
   int pos = 0, pos_prev =0;
   int len_max = 0, len_new = 0, len_cum = 0;
   int i;
   char Lmer[seq-> L];
   memcpy(Lmer,seq->line2,seq -> L);
   Lmer_sLmer(Lmer,seq->L);
   for (i = 0 ; i < seq -> L; i++){
      len_cum++;
      if(Lmer[i] > '\003' ){
         pos_prev = i + 1 - len_cum;
         len_new = i - pos_prev;
         len_cum = 0; 
         if(len_new > len_max){
            pos = pos_prev;
            len_max = len_new;
         }
      }
   }
   if (len_cum == seq -> L){
      return 1;
   }
   // Consider last piece
   if(len_cum > len_max){
      len_max = len_cum ;
      pos = seq -> L  - len_max;
   }
   if( len_max < minL){
      return 0;
   }else{
      memmove(seq -> line2 ,seq -> line2+pos,len_max);
      seq -> line2[len_max] = '\0';
      memmove(seq -> line4 ,seq -> line4+pos,len_max);
      seq -> line4[len_max] = '\0';
      seq -> L = len_max;
      char add[10]; 
      sprintf(add," TRIMN:%d:%d",pos,pos+len_max);
      strcat(seq -> line3, add);
      return 2; 
   }

}

/* returns 1 if no N's found, 2 if trimmed*/
int Ntrim_ends(Fq_read *seq, int minL){
   char Lmer[seq-> L];
   int t_start = 0 ;  
   int t_end = seq -> L - 1;
   memcpy(Lmer,seq->line2,seq -> L);
   Lmer_sLmer(Lmer,seq->L);
   while (Lmer[t_start] > '\003')
      t_start++;
   while (Lmer[t_end] > '\003')
      t_end--;
   if( (t_end - t_start) == (seq -> L - 1) )
      return 1; 
   else if ((t_end - t_start) < minL -1 ){
      return 0; 
   } else{
     (seq -> L) = t_end - t_start + 1;
     memmove(seq -> line4, seq -> line4 + t_start, seq -> L);
     memmove(seq -> line2, seq -> line2 + t_start, seq -> L);
     seq -> line4[seq -> L] = '\0';
     seq -> line2[seq -> L] = '\0';
     char add[10]; 
     sprintf(add," TRIMN:%d:%d",t_start, t_end);
     strcat(seq -> line3, add);
     return 2;
   }
} 



/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceQ(Fq_read *seq, int minQ, bool trimQ, int minL){
   // Beginning of sequence:
   int L = (seq->L)-1;
   int t_start = 0;
   int t_end = L;
   int i;
   while ( seq->line4[t_start] < (ZEROQ + minQ) )
      t_start++;
   while ( seq->line4[t_end] < (ZEROQ + minQ) )
      t_end--;
   for ( i = t_start ; i <  t_end ; i++){
      if ( seq -> line4[i] < (ZEROQ + minQ)){
         return 0 ; 
      }
   }
   if ((t_end - t_start) == L ){
      // accepting sequence as is
      return 1; 
   } 
   if(!trimQ){
      // if !trimQ, discard the sequence if not complete
      return 0;
   }
   // trim the sequence
   if ((t_end - t_start) < minL -1  ){
      // Ignoring sequence
      return 0; 
   } 
   (seq -> L) = t_end - t_start + 1;
   memmove(seq -> line4, seq -> line4 + t_start, seq -> L);
   memmove(seq -> line2, seq -> line2 + t_start, seq -> L);
   seq -> line4[seq -> L] = '\0';
   seq -> line2[seq -> L] = '\0';
   char add[10]; 
   sprintf(add," TRIMQ:%d:%d",t_start, t_end);
   strcat(seq -> line3, add);
   return 2;
}
/* Checking N  base callings*/
/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceN(Fq_read *seq, int trimN, int minL ){
   return (trimN == 0)? no_N(seq)  : 
          (trimN == 1)? Nfree_Lmer(seq, minL):
          (trimN == 2)? Ntrim_ends(seq, minL): -1;
}

// Check if a read is in the sequence or its reverse complement
bool is_read_in_seq(Node *tree, Fq_read* seq, int L){
   char Lmer[seq -> L];
   memcpy(Lmer,seq->line2,seq -> L);
   // translate to 0,1,2,3
   Lmer_sLmer(Lmer,seq -> L);
   // check path
   if (check_path(tree, Lmer, L, seq -> L)){
         return true;
   } else {
        rev_comp(Lmer,seq -> L);
        return (check_path(tree, Lmer, L, seq -> L));
   }
   
}
