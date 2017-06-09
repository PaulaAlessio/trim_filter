#include "bloom_maker.h"
#include "tree_fasta.h"
#include <string.h>
#include <stdio.h>


// Global variables (lookup table)
// Used to compactify kmers
static uint8_t fw0[256], fw1[256], fw2[256], fw3[256]; 
static uint8_t bw0[256], bw1[256], bw2[256], bw3[256]; 

// Initialize Lookup tables
void init_LUTs(){

   memset(fw0,0xFF,256);
   memset(fw1,0xFF,256);
   memset(fw2,0xFF,256);
   memset(fw3,0xFF,256);
   memset(bw0,0xFF,256);
   memset(bw1,0xFF,256);
   memset(bw2,0xFF,256);
   memset(bw3,0xFF,256);
   // Mappings 
   fw0['a'] = 0x00; fw0['c'] = 0x40; fw0['g'] = 0x80; fw0['t'] = 0xC0; 
   fw0['A'] = 0x00; fw0['C'] = 0x40; fw0['G'] = 0x80; fw0['T'] = 0xC0;
   fw1['a'] = 0x00; fw1['c'] = 0x10; fw1['g'] = 0x20; fw1['t'] = 0x30; 
   fw1['A'] = 0x00; fw1['C'] = 0x10; fw1['G'] = 0x20; fw1['T'] = 0x30;
   fw2['a'] = 0x00; fw2['c'] = 0x04; fw2['g'] = 0x08; fw2['t'] = 0x0C; 
   fw2['A'] = 0x00; fw2['C'] = 0x04; fw2['G'] = 0x08; fw2['T'] = 0x0C;
   fw3['a'] = 0x00; fw3['c'] = 0x01; fw3['g'] = 0x02; fw3['t'] = 0x03; 
   fw3['A'] = 0x00; fw3['C'] = 0x01; fw3['G'] = 0x02; fw3['T'] = 0x03;

   bw0['a'] = 0xC0; bw0['c'] = 0x80; bw0['g'] = 0x40; bw0['t'] = 0x00;
   bw0['A'] = 0xC0; bw0['C'] = 0x80; bw0['G'] = 0x40; bw0['T'] = 0x00;
   bw1['a'] = 0x30; bw1['c'] = 0x20; bw1['g'] = 0x10; bw1['t'] = 0x00;
   bw1['A'] = 0x30; bw1['C'] = 0x20; bw1['G'] = 0x10; bw1['T'] = 0x00;
   bw2['a'] = 0x0C; bw2['c'] = 0x08; bw2['g'] = 0x04; bw2['t'] = 0x00;
   bw2['A'] = 0x0C; bw2['C'] = 0x08; bw2['G'] = 0x04; bw2['T'] = 0x00;
   bw3['a'] = 0x03; bw3['c'] = 0x02; bw3['g'] = 0x01; bw3['t'] = 0x00;
   bw3['A'] = 0x03; bw3['C'] = 0x02; bw3['G'] = 0x01; bw3['T'] = 0x00;

}


void init_Bfilter(int kmersize, uint64_t bfsizeBits, int hashNum, Bfilter *bf){
   bf -> kmersize = kmersize;
   bf -> kmersizeBytes = (kmersize + 4 - 1) / 4;
   bf -> hashNum = hashNum;
   bf -> bfsizeBits = bfsizeBits;
   bf -> bfsizeBytes = bfsizeBits/BITSPERCHAR;
   bf -> filter = (unsigned char *) malloc (bf ->  bfsizeBytes * sizeof(unsigned char));
   memset(bf -> filter, 0, bf -> bfsizeBytes);
}

void free_Bfilter(Bfilter * bf){
   free(bf -> filter);
}

void init_procs(Procs_kmer *procs, int kmersize, int hashNum){

   procs -> kmersize = kmersize; 
   procs -> kmersizeBytes = kmersize / 4; 
   procs -> halfsizeBytes = kmersize / 8; 
   procs -> hangingBases = 0;
   procs -> hasOverhead = 0;
   procs -> hashNum = hashNum;
   if (kmersize % 8 != 0) {
         procs -> halfsizeBytes++; 
         if( (procs -> hangingBases = kmersize % 4) > 0 ){
               procs -> kmersizeBytes++;
               procs -> hasOverhead = 1; 
         }
     }

     procs -> compact = (unsigned char*)malloc(procs -> kmersizeBytes*sizeof(unsigned char));
     procs -> hashValues = (uint64_t*)malloc(hashNum*sizeof(uint64_t));

}

void free_procs(Procs_kmer *procs){
   free(procs->compact);
   free(procs->hashValues);
}

static int compact_kmer(const unsigned char *sequence, uint64_t position, Procs_kmer *procs){
   unsigned char * m_fw, *m_bw;
   m_fw = (unsigned char*)malloc(procs -> kmersizeBytes*sizeof(unsigned char));
   m_bw = (unsigned char*)malloc(procs -> kmersizeBytes*sizeof(unsigned char));
   memset(m_fw, 0, procs -> kmersizeBytes);
   memset(m_bw, 0, procs -> kmersizeBytes);
   int idx;   // indexes the compactified 
   uint64_t b = position, revb = position + procs->kmersize - 1; // position in the sequence
   for (idx = 0 ; idx < procs -> halfsizeBytes; idx++){
        m_fw[idx] |= fw0[sequence[b++]];
        m_fw[idx] |= fw1[sequence[b++]];
        m_fw[idx] |= fw2[sequence[b++]];
        if ( (m_fw[idx] == 0xFF) || ( fw3[sequence[b]] == 0xFF ) ){
            // discards kmers with N's
            free(m_fw); free(m_bw);
            return 0;
        }
        m_fw[idx] |= fw3[sequence[b++]];
        
        m_bw[idx] |= bw0[sequence[revb--]];
        m_bw[idx] |= bw1[sequence[revb--]];
        m_bw[idx] |= bw2[sequence[revb--]];
        if ( (m_bw[idx] == 0xFF) || ( bw3[sequence[revb]] == 0xFF ) ){
            // discards kmers with N's
            free(m_fw); free(m_bw);
            return 0;
        }

        m_bw[idx] |= bw3[sequence[revb--]];
        // Check which one is lexicographically smaller 
        if(m_fw[idx] < m_bw[idx] ){ // go on with forward
           free(m_bw);
           for(++idx; idx < (procs-> kmersizeBytes - procs -> hasOverhead) ; idx++ ){
                m_fw[idx] |= fw0[sequence[b++]];
                m_fw[idx] |= fw1[sequence[b++]];
                m_fw[idx] |= fw2[sequence[b++]];
                if ( (m_fw[idx] == 0xFF) || ( fw3[sequence[b]] == 0xFF ) ){
                    // discards kmers with N's
                    free(m_fw); 
                    return 0;
                }
                m_fw[idx] |= fw3[sequence[b++]];
           }
           if(procs -> hasOverhead){
              idx++;
              switch (procs->hangingBases){
                case 1: 
                  m_fw[idx] |= fw0[sequence[b++]];
                  break;
                case 2: 
                  m_fw[idx] |= fw0[sequence[b++]];
                  m_fw[idx] |= fw1[sequence[b++]];
                  break;
                case 3: 
                  m_fw[idx] |= fw0[sequence[b++]];
                  m_fw[idx] |= fw1[sequence[b++]];
                  m_fw[idx] |= fw2[sequence[b++]];
                  break;
                default: 
                   printf("ERROR\n");
                   return 0; 
              }
              if ( (m_fw[idx] == 0xFF)  ){
                    // discards kmers with N's
                    free(m_fw); 
                    return 0;
              }
           }
           memcpy(procs -> compact, m_fw, procs->kmersizeBytes);  
           free(m_fw);
           return 1;  
        } else if ( m_fw[idx] > m_bw[idx] ){ // go on with backward
           free(m_fw);
           for(++idx; idx < (procs-> kmersizeBytes - procs -> hasOverhead) ; idx++ ){
                m_bw[idx] |= bw0[sequence[revb--]];
                m_bw[idx] |= bw1[sequence[revb--]];
                m_bw[idx] |= bw2[sequence[revb--]];
                if ( (m_bw[idx] == 0xFF) || ( bw3[sequence[revb]] == 0xFF ) ){
                    // discards kmers with N's
                    free(m_bw); 
                    return 0;
                }
                m_bw[idx] |= bw3[sequence[revb--]];
           }
           if(procs -> hasOverhead){
              idx++;
              switch (procs->hangingBases){
                case 1: 
                  m_bw[idx] |= bw0[sequence[revb--]];
                  break;
                case 2: 
                  m_bw[idx] |= bw0[sequence[revb--]];
                  m_bw[idx] |= bw1[sequence[revb--]];
                  break;
                case 3: 
                  m_bw[idx] |= bw0[sequence[revb--]];
                  m_bw[idx] |= bw1[sequence[revb--]];
                  m_bw[idx] |= bw2[sequence[revb--]];
                  break;
                default: 
                   printf("ERROR\n");
                   return 0; 
              }
              if ( (m_bw[idx] == 0xFF)  ){
                    // discards kmers with N's
                    free(m_bw); 
                    return 0;
              }
           }

           memcpy(procs -> compact, m_bw, procs->kmersizeBytes);  
           free(m_bw);
           return 2;
        }
        
   }
   // If it is a palindrome... We have to cover this case!
   free(m_bw);
   for(; idx < (procs -> kmersizeBytes - procs -> hasOverhead); idx++ ){
        m_fw[idx] |= fw0[sequence[b++]];
        m_fw[idx] |= fw1[sequence[b++]];
        m_fw[idx] |= fw2[sequence[b++]];
        if ( (m_fw[idx] == 0xFF) || ( fw3[sequence[b]] == 0xFF ) ){
            // discards kmers with N's
            free(m_fw); free(m_bw);
            return 0;
        }
        m_fw[idx] |= fw3[sequence[b++]];
   } 
   // finish last byte     
   if(procs -> hasOverhead){
   idx++;
   switch (procs->hangingBases){
     case 1: 
       m_fw[idx] |= fw0[sequence[b++]];
       break;
     case 2: 
       m_fw[idx] |= fw0[sequence[b++]];
       m_fw[idx] |= fw1[sequence[b++]];
       break;
     case 3: 
       m_fw[idx] |= fw0[sequence[b++]];
       m_fw[idx] |= fw1[sequence[b++]];
       m_fw[idx] |= fw2[sequence[b++]];
       break;
     default: 
        printf("ERROR\n");
        return 0; 
   }
      if ( (m_fw[idx] == 0xFF)  ){
         // discards kmers with N's
         free(m_fw); 
         return 0;
      }
   }
   memcpy(procs -> compact, m_fw, procs->kmersizeBytes);  
   free(m_fw);
   return 3;
   
}


static void multiHash( Procs_kmer* procs){
   int kmerSizeInBytes = ((procs -> kmersize + 4 - 1)/4);
   int i; 
   for(i = 0; i <  procs -> hashNum; i++){
      procs -> hashValues[i] = CityHash64WithSeed( (const char *) procs->compact,
          kmerSizeInBytes, i );
   }
}

static int insert_and_fetch(Bfilter *bf, Procs_kmer* procs){
   int result = 1;
   int i = 0;
   uint64_t modValue; 
   //iterates through hashed values adding it to the filter
   for (i = 0; i < procs -> hashNum; i++) {
      modValue = (procs -> hashValues[i]) % (bf -> bfsizeBits);
//      printf("%ld %ld\n", procs -> hashValues[i] ,modValue);
      result &= (__sync_fetch_and_or(  &(bf -> filter)[modValue/BITSPERCHAR], 
              bitMask[modValue % BITSPERCHAR] ) ) >> (modValue % BITSPERCHAR) & 1; 
   }
   return result;
}

static int  contains(Bfilter *bf, Procs_kmer* procs){
   int i = 0;
   uint64_t modValue; 
   //iterates through hashed values and check whether they are in the filter
   for (i = 0; i < procs -> hashNum; i++) {
      modValue = (procs -> hashValues[i]) % (bf -> bfsizeBits);
      unsigned char bit = bitMask[modValue % BITSPERCHAR];
      if (( (bf -> filter)[modValue / BITSPERCHAR] & bit) != bit) {
          return 0;
      }
   }
   return 1; 
   
}


double is_read_in_filter(unsigned char *read, int L, Procs_kmer *ptr_procs, Bfilter *bf){
   int position; 
   int maxN = L - bf -> kmersize;
   double score = 0;
   int stretch = 0; 
   for (position = 0; position <  maxN; position++ ){
      
     if( compact_kmer(read,position,ptr_procs)) {
         multiHash( ptr_procs);
         if(contains( bf, ptr_procs)){
            stretch++;
            score += 1 - 1.0/(2*stretch);
         } else {
            stretch = 0;       
         }
    
     } else{
       stretch = 0; 
     }
   }
   return score; 
}


int create_bloom_filter(char **fastafiles, int Nfiles, Bfilter *bf ){
   int i_files, i_entries;
   int position;
   init_LUTs();
   Procs_kmer procs;
   init_procs(&procs,bf -> kmersize, bf->hashNum);
   int isvalid; 
   for (i_files = 0 ; i_files < Nfiles; i_files++ ){
      Fasta *fasta_seq; 
      int Nentries = read_fasta(fastafiles[i_files], &fasta_seq);
      for (i_entries = 0; i_entries < Nentries; i_entries++){
         int maxN = fasta_seq[i_entries].N - bf->kmersize;
         for (position = 0 ; position < maxN; position++){
            isvalid = compact_kmer(fasta_seq[i_entries].seq, position, &procs);
            if(isvalid > 0){
               multiHash( &procs);
               insert_and_fetch( bf, &procs);
            //   int result = insert_and_fetch( bf, &procs);
            }
//           int k, si ; 
//            if (position < 10 ){
//              for (si = 0 ; si  < sizeof(procs.compact); si++){
//                 for(k = 7 ; k >=0; k--){
//                    printf("%d", (procs.compact[si] & (1<< k)) >> k );
//                 } 
//              }
//               printf("\n");
//            
//            }
         }
      }   
      free_fasta(&fasta_seq,Nentries);
   }
   free_procs(&procs);
   return 1; 
}


void save_bloomfilter(char *outputfile, Bfilter *bf){
   fprintf(stderr, "Store the bloom filter in: %s\n", outputfile);
   fprintf(stderr, "Bloom filter size in bytes, %ld\n", bf->bfsizeBytes);
   FILE *fout = fopen(outputfile, "wb");
   fwrite(bf->filter,sizeof(char), bf -> bfsizeBytes, fout);
   fclose(fout);
}
