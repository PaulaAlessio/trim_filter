//defines.h
/* contains preprocessor directives*/
#ifndef DEFINES_H
#define DEFINES_H

#define FASTA_ENTRIES 1
#define B_LEN 131072
#define ACGT 4
#define BUF_NODE 2000
#define POOL_NUM 10 
#define PATH_MAXL 300


#define bool int
#define true 1
#define false 0

#ifndef max
   #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
   #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef mem_usage
   #define mem_usageMB()  fprintf(stderr,"- Current allocated memory: %ld MB.\n",alloc_mem >> 20)
#endif

#ifndef mem_usage
   #define mem_usage()  fprintf(stderr,"- Current allocated memory: %ld Bytes.\n",alloc_mem)
#endif

#endif

