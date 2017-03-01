#include "init.h"
#include <string.h> /* for strcmp*/
#include <stdio.h> /* for printf*/
#include <stdlib.h> /* for exit*/
#include <getopt.h> /* for getopt*/

void printHelpDialog(){
   const char dialog[] =
     "Usage: trimFilter --ifq INPUT.fq -l READLEN [-q MINQ | --ifa INPUT.fa | --idx IDX_FILE |" 
                   " -o O_PREFIX | (--tree|-t) | --trimQ | --trimN [all|ends] | (--help|-h) ]\n"
     "Reads in a fq file (gz, bz2, z formats also accepted) and removes: \n"
     "  * Low quality reads (trim reads passing the --trim option),\n"
     "  * Reads containing N base callings,\n"
     "Outputs 4 [O_PREFIX]_fq.gz files containing: \"good\" reads, discarded low Q reads,\n"
     "discarded reads containing N's, discarded contaminations.\n"
     "Options:\n"
     " -h, --help   prints help dialog.\n"
     " -l           Read length: length of the reads. Required option.\n"
     " -f, --ifq    Fastq input file [*fq|*fq.gz|*fq.bz2]: Required option.\n"
     " -a, --ifa    Fasta input file [*fa|*fa.gz|*fa.bz2]: Optional.\n"
     " -x, --idx    Index input file: Optional.\n"
     " -k, --lmer   Lmer length: Length of the Lmers to be searched on the contaminations *fa.\n"
     "              Default: length of the reads.\n"
     " -t, --tree   Construct a tree for the search process. If missing, SA constructed.\n"
     " -Q, --trimQ  no (or no flag): does nothing to low quality reads,\n"
     "              ends: trims the ends of the read if their quality is below the threshold -q,\n"
     "              all:  removes all reads containing at least one low quality nucleotide.\n"
     " -N, --trimN  no (or no flag): does nothing to reads containing N's,\n" 
     "              all: removes all reads containing N's,\n"
     "              ends: trims ends of reads with N's,\n"
     "              strips: looks for the largest substring containing no N's.\n"
     " -q           Minimum quality allowed. Optional option. Default 27 .\n"
     " -o           Output prefix (containing the path). Optional option. Default ./out.\n";
   fprintf(stderr, "%s", dialog );
   exit(EXIT_SUCCESS);
}

static Param default_par(){
    Param par; 
    par.Ifq = "";    
    par.Ifa = "";    
    par.Iidx = "";    
    par.Oprefix = "out";
    par.trimQ = 0;
    par.trimN = 0;
    par.is_fa = false; 
    par.is_idx = false; 
    par.tree = false; 
    par.minQ = 27; 
    par.L = 0;  
    par.Lmer_len = 0;  
    return par; 
}
 
Param  get_arg( int argc, char **argv){
   static struct option long_options[] = {
      {"ifq",required_argument,0,'f'},  
      {"ifa",required_argument,0,'a'},  
      {"idx",required_argument,0,'x'},  
      {"lmer",required_argument,0,'k'},  
      {"trimQ",required_argument,0,'Q'},  
      {"trimN",required_argument,0,'N'},  
      {"tree",no_argument,0,'t'}, 
      {"help",no_argument,0,'h'}
   };
   int tmp;
   Param par = default_par(); 
   while ( (tmp = getopt_long(argc, argv,"hl:q:o:tf:a:x:k:Q:N:", long_options, 0)) != -1){
      switch (tmp){
         case 'h': 
            printHelpDialog();
            break;
         case 'l':
            par.L = atoi(optarg); 
            break;    
         case 'q':
            par.minQ = atoi(optarg);
            break;
         case 'o':
            par.Oprefix = optarg; 
            break;
         case 't': 
            par.tree = true;
            break;
         case 'f':
            par.Ifq =  optarg;
            break;
         case 'a':
            par.Ifa =  optarg;
            par.is_fa = true;
            break;
         case 'x':
            par.Iidx = optarg;
            par.is_idx = true;
            break;
         case 'k':
            par.Lmer_len = atoi(optarg);
            break;
         case 'Q':
            if (!strcmp("no",optarg))
               par.trimQ = NO;
            else if(!strcmp("all",optarg))
               par.trimQ = ALL;
            else if (!strcmp("ends",optarg))
               par.trimQ = ENDS;
            else{
               fprintf(stderr, "OPTIONS_ERROR: trimQ accepts: no, all, end as arguments.\n");
               printHelpDialog();
            }

            break;
         case 'N':
            if (!strcmp("no",optarg))
               par.trimN = NO;
            else if(!strcmp("all",optarg))
               par.trimN = ALL;
            else if (!strcmp("ends",optarg))
               par.trimN = ENDS;
            else if (!strcmp("strip",optarg))
               par.trimN = STRIP;
            else{
               fprintf(stderr, "OPTIONS_ERROR: trimN accepts: no, all, end, strip as arguments.\n");
               printHelpDialog();
            }
            break;
         default:
            fprintf(stderr,"?? getopt returned character code 0%o ??\n", tmp);
            printHelpDialog();
            break;
      
      }
   }
   // Check mandatory options: 
   if( !strcmp(par.Ifq,"") || (par.L==0) ){
        fprintf(stderr,"OPTIONS_ERROR: --ifq and -l are mandatory options.\n");
        printHelpDialog();
   }
   // Check not yet implemented functions: 
   if( !par.tree ){
        fprintf(stderr,"OPTIONS_WARNING: SA option has not been implemented. Pass --tree.\n" );
        printHelpDialog();
   }
   fprintf(stderr, "- Fastq input file: %s \n",par.Ifq);
   fprintf(stderr, "- Read length: %d.\n",par.L);
   fprintf(stderr, "- Min Quality: %d%s.\n",par.minQ,(par.minQ==27)?" (default)":"");
   fprintf(stderr, "- Output prefix: %s%s.\n",par.Oprefix,
                      strcmp(par.Oprefix,"out")?"":" (default)");
   (par.Lmer_len > 0)?
      fprintf(stderr,"- Verifying matches using Lmers of length %d.\n", par.Lmer_len):
      fprintf(stderr,"- No Lmer_len specified. Using the length of the sequence L = %d.\n", 
            par.Lmer_len = par.L);
   fprintf(stderr, "- Triming N bases: %s.\n", (par.trimN == NO )? "no": 
               (par.trimN == ALL)? "removing all reads with N" : 
               (par.trimN == ENDS)? "trimming N's at the ends" : 
               "stripping reads, extracting the longest N free subsequence");
   fprintf(stderr, "- Triming low Q bases: %s.\n", (par.trimQ == NO) ? "no":
               (par.trimQ == ALL)? "removing all reads with low Q nucleotides":
               "trimming low Q nucleotides at the ends");
   // Checking for fasta file
   par.is_fa?
           fprintf(stderr, "- Removing contaminations in %s.\n", par.Ifa): 
           fprintf(stderr, "- NOT Removing contaminations (No fasta file passed).\n");
   // Check the tree option (SA not implemented yet!)
   par.tree?
      fprintf(stderr, "- Performing a quaternary tree search.\n"): 
      fprintf(stderr, "- Using SA and a binary search.\n");
   // Check whethere an index file was available (Option not implemented yet!)
   if(par.is_idx){
      fprintf(stderr, "- Reading %s from file %s\n",par.tree?"tree":"SA",par.Iidx);
      fprintf(stderr,"OPTIONS_WARNING: tree file option has not been implemented yet.\n");
      printHelpDialog();
   }else{
      fprintf(stderr, "- No %s file passed. Computing it on the fly.\n",
              par.tree?"tree":"SA");
   }
   return par; 
}
  
