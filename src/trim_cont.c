#include <string.h>
#include "trim_cont.h"
#include "Lmer.h"
#include "defines.h"

/* returns 0 if no N's found, 1 if N's found*/
int no_N(Fq_read *seq){
   int i;
   char Lmer[seq-> L];
   memcpy(Lmer,seq->line2,seq -> L);
   Lmer_sLmer(Lmer,seq->L);
   for(i = 0; i < seq->L; i++){
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

/* 0 if sequence contains low quality nucleotides, 1 otherwise*/
int no_lowQ(Fq_read *seq, int minQ) {
   int i;
   for (i = 0; i < seq -> L; i++ ){
     if (seq->line4[i] < (ZEROQ + minQ))
         return 0;  
   }
   return 1;
}

/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int Qtrim_ends(Fq_read *seq, int minQ, int minL){
   // Beginning of sequence:
   int L = (seq->L)-1;
   int t_start = 0;
   int t_end = L;
   while ( seq->line4[t_start] < (ZEROQ + minQ) )
      t_start++;
   while ( seq->line4[t_end] < (ZEROQ + minQ) )
      t_end--;
   if ((t_end - t_start) == L ){
      // accepting sequence as is
      return 1; 
   } 

   // trim the sequence
   if ((t_end - t_start) < minL - 1 ){
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

/* 0 if not used, 1 if accepted as is*/
int Qtrim_frac(Fq_read *seq, int minQ, int nlowQ ){
   int ilowQ = 0, i = 0; 
   while (i < (seq->L) &&  ilowQ < nlowQ){
      if(seq->line4[i] < (ZEROQ + minQ)) ilowQ++;
      i++;
   }  
   return (ilowQ == nlowQ) ? 0: 1; 
}

/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int Qtrim_endsfrac(Fq_read *seq, int minQ, int minL, int nlowQ ){
   // Beginning of sequence:
   int L = (seq->L)-1;
   int t_start = 0;
   int t_end = L;
   int ilowQ = 0;
   int i;
   while ( seq->line4[t_start] < (ZEROQ + minQ) )
      t_start++;
   while ( seq->line4[t_end] < (ZEROQ + minQ) )
      t_end--;
   i = t_start; 
   while (i < t_end && ilowQ < nlowQ){
      if(seq->line4[i] < (ZEROQ + minQ)) ilowQ++;
      i++;
   }  
   if(ilowQ == nlowQ ) return 0; 
   if ((t_end - t_start) == L ){
      // accepting sequence as is
      return 1; 
   } 

   // trim the sequence
   if ((t_end - t_start) < minL - 1 ){
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

/* returns 2, since they are all accepted and trimmed*/
int Qtrim_global(Fq_read *seq, int left, int right ){
   int t_end = seq ->L - left;
   (seq -> L) -= (left+right);
   memmove(seq -> line4, seq -> line4 + left, seq -> L);
   memmove(seq -> line2, seq -> line2 + left, seq -> L);
   seq -> line4[seq -> L] = '\0';
   seq -> line2[seq -> L] = '\0';
   char add[10]; 
   sprintf(add," TRIMQ:%d:%d",left, t_end);
   strcat(seq -> line3, add);
   return 2;
   
}

/* Checking N  base callings*/
/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceN(Fq_read *seq, Param par ){
   return (par.trimN == NO)? 1 : 
          (par.trimN == ALL)? no_N(seq):
          (par.trimN == ENDS)? Ntrim_ends(seq, par.minL): 
          (par.trimN == STRIP)? Nfree_Lmer(seq, par.minL): -1;
}



/* 0 if not used, 1 if accepted as is, 2 if accepted and trimmed*/
int trim_sequenceQ(Fq_read *seq, Param par){
   return (par.trimQ == NO)? 1 :
          (par.trimQ == ALL)? no_lowQ(seq, par.minQ): 
          (par.trimQ == ENDS) ? Qtrim_ends(seq,par.minQ, par.minL): 
          (par.trimQ == FRAC) ? Qtrim_frac(seq,par.minQ, par.nlowQ): 
          (par.trimQ == ENDSFRAC) ? Qtrim_endsfrac(seq, par.minQ, par.minL, par.nlowQ): 
          (par.trimQ == GLOBAL) ? Qtrim_global(seq,par.globleft,par.globright): -1;  
          
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
