// fopen_gen.h 
/* This functions allow us to open different types of 
 * compressed files: (*.tar.gz, *.gz, *.bz2, *.bam, etc)
 * using the function fopen_gen. It returns a pointer 
 * to FILE that can be then used as usual with fread. */
#ifndef fopen_gen_H
#define fopen_gen_H

#include <stdio.h>

#define READ_END 0
#define WRITE_END 1
#define PERMISSIONS 0640


#ifdef __STDC__
FILE* fdopen( int, const char* );
#endif

int setCloexec(int fd);
FILE* fopen_gen(const char *path, const  char * mode);

/* Static functions
 * fopen_gen in READ mode:
 *  static const char* zcatExec(const char* path);
 *  static int uncompress(const char* path);
 *  static FILE* funcompress(const char* path);
 * fopen_gen in WRITE mode:
 * static const char* catExec(const char* path);
 * static int compress(const char* path);
 * static FILE* compress(const char* path);
*/

#endif
