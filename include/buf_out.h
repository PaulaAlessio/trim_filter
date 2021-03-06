// buf_out.h 
/* function that buffers the output 
 * before writing to disk
 */
#ifndef BUF_OUT_H
#define BUF_OUT_H

#define GOOD 0
#define CONTAMINATION 1
#define LOWQ 2
#define NNNN 3

void buffer_output(FILE *fout, const char *a, const int len, const int fd_i);

void write_summary(int nreads, int n_good, int n_NNNN, int n_lowq, 
      int n_trimN, int n_trimQ, int n_cont, char *summary);
#endif
