# trimFilter user manual

 This programs reads a `fastq` file as an input, and performs the following filters:
 - If a `fasta` file is passed as an input (e.g., rRNA ), it will discard reads
   that are contained in it. 
 - Discards/trims low quality reads (see below for a description of  *low quality*).
 - Discards/trims reads containing N base calls. (se dice as'i?)

## Installation

 - Requires: fill that in. Basically: Linux, gcc  
 - Clone the repository.
 - Run make.
  
## Running the program

`C` executable (in folder `bin`):

```
Usage: trimFilter --ifq INPUT.fq -l READLEN [-q MINQ | --ifa INPUT.fa | --idx IDX_FILE |
                  -o O_PREFIX | (--tree|-t) | --trimQ [no|all|ends|frac|endsfrac|global] | 
                  (--percent percent) | (--global n1 n2) | --trimN [all|ends] | (--help|-h) ]
Reads in a fq file (gz, bz2, z formats also accepted) and removes: 
  * Low quality reads (trim reads passing the --trim option),
  * Reads containing N base callings,
Outputs 4 [O_PREFIX]_fq.gz files containing: "good" reads, discarded low Q reads,
discarded reads containing N's, discarded contaminations.
Options:
 -h, --help    prints help dialog.
 -l            Read length: length of the reads. Required option.
 -q            Minimum quality allowed (int). Optional option. Default 27 .
 -o            Output prefix (containing the path). Optional option. Default ./out.
 -f, --ifq     Fastq input file [*fq|*fq.gz|*fq.bz2]: Required option.
 -a, --ifa     Fasta input file [*fa|*fa.gz|*fa.bz2]: Optional.
 -x, --idx     Index input file: Optional.
 -k, --lmer    Lmer length: Length of the Lmers to be searched on the contaminations *fa.
               Default: length of the reads.
 -t, --tree    Construct a tree for the search process. If missing, SA constructed. Optional
 -Q, --trimQ   no (or no flag): does nothing to low quality reads,
               all:  removes all reads containing at least one low quality nucleotide.
               ends: trims the ends of the read if their quality is below the threshold -q,
               frac: discards a read if the fraction of bases whose quality  
                     lies below the threshold -q is over 5% or a user defined 
                     percentage in -p.
               endsfrac: trims the ends and then discards the read if there are more
               global: removes n1 cycles on the left and n2 on the right, user 
                       specified in -g.
 -p, --percent percentage of low quality bases to be admitted before discarding a read, 
 -g, --global  required option if --trimQ global is passed. two integers, 
               n1, n2, have to be passed specifying the number of cycles 
               to be globally cut from the left and right, respectively.
 -N, --trimN   no (or no flag): does nothing to reads containing N's,
               all: removes all reads containing N's,
               ends: trims ends of reads with N's,
               strips: looks for the largest substring containing no N's.
```

`Rmd` script (in folder `R`):

```
Usage: 

Rscript -e "rmarkdown::render('PATH/TO/summary_report.Rmd', 
            params=list(inputfolder='O_PREFIX/FOLDER/'),
            output_file='PATH/TO/HTML_OUTPUT_FILE.html')"
```


## Output description

- `O_PREFIX_good.fq.gz`: contains reads that passed all filters.
- `O_PREFIX_NNNN.fq.gz`: contains reads with *too many* N's. 
- `O_PREFIX_lowQ.fq.gz`: contains reads discarded due to low quality issues.
- `O_PREFIX_cont.fq.gz`: contains contamination reads (matching the fasta sequence). 
- `O_PREFIX_summary.bin`: binary file containing 7 `int`'s: 
    * # initial reads, 
    * # accepted reads, 
    * # reads with N's that were discarded, 
    * # low Quality reads that were discarded, 
    * # reads trimmed due to the presence of N's, 
    * # reads trimmed due to low Quality, 
    * # reads identified as contaminations. 


## Filters

#### Impurities  

 Contaminations are removed if a fasta file is given as an input. 
 Two methods have been implemented to check for contaminations: 

- **tree**: this method is designed to identify impurities from 
  small sequences, such as 
  if `--tree` is passed, then a tree with all 
  possible substrings of length `Lmer` is constructed,
  set to the read length by default. If all `Lmer`'s 
  of a read are found in the tree, then the read is 
  classified as a contamination and redirected to 
  `*_cont.fq.gz`. Some considerations about this 
  method: 
    - it does not introduce false negatives, i.e., it cannot
      be that a read is exactly contained in the `*fasta` file
      and is not discarded. 
    - is very fast. Every search is `O( 2 * Lmer * (L - Lmer + 1))` (the reverse complement
      of the read has to be checked also).
    - has the drawback that is very memory intensive.  

- **SA**: if the `fasta` file against which we have to 
  look for contaminations is large (e.g. the *Drosophila* genome), 
  then constructing a tree is not viable due to memory limitations. 
  In this case, we construct a suffix array from the genome, that 
  should not exceed 4gb pairs. The strategy followed is analogous
  to the one in `STAR`: a binary search for exact matches is then
  performed, but a prefix index of `m` bases is constructed, 
  so that we do not have to jump so much within a binary search. 
  Some remarks: 
   - it is not fully implemented yet. It works for a 
     `fasta` file with one entry. 
   - The complexity is `O(log(Lgenome)*2*L)`, where L is the read length. 
 

#### LowQ

- `--trimQ no` or flag absent: nothing is done to the reads with low quality.
- `--trimQ all`: all reads containing at least one low quality nucleotide are
  redirected to  `*_lowq.fq.gz`
- `--trimQ ends`: look for low quality (below minQ) base callings at the beginning and at the end of the read.
  Trim them at both ends until the quality is above the threshold. Keep the read in `*_good.fq.gz`
  and annotate in the fourth line where the read has been trimmed (starting to count
  from 0) if the length of the remaining part is larger than the half of the original read length.
  Redirect the read to `*_lowq.fq.gz` otherwise.
    Examples (-q 27 [<]): 
 ```
 @ read 1081133
 GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA  GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGG
 +position: -144740                                  +position: -144740 TRIMQ:0:48
 IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 3435                                         @ read  3435 
 CAGTTCTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTGACCCA  CTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTG
 +position: -19377                                   +position: -19377  TRIMQ:5:44
 99999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99999  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 108110                                       @ read 108110 
 CAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAAC  AGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAA
 +position: -4142524                                 +position: -4142524  TRIMQ:1:48
 9IIIIIIIIIIIIIIIIIII9IIII9III9IIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIIIIIIIII9IIII9III9IIIIIIIIIIIIIIIIIII
 @ read 108111                                      
 AGCGTGACGGTGTAACGCCCGCCTTTTGATGACTGGGTTTCAAAGAAACG  
 +position: -3336785                                 
 999999999999999IIIIIIII9IIIIIIIII9II99999999999999  
 ```

- `--trim frac (--percent p)`: redirect  the read to `*_lowq.fq.gz` if there are more
   than `p%` nucleotides whose quality lies below the threshold. 
   `p=5` per default.  
    Examples (-q 27 [<]): 
 ```
 @ read 1081133
 GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA  GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA
 +position: -144740                                  +position: -144740 
 IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9
 @ read 3435                                         
 CAGTTCTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTGACCCA  
 +position: -19377                                   
 99999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99999  
 @ read 108110                                       
 CAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAAC  
 +position: -4142524                                 
 9IIIIIIIIIIIIIIIIIII9IIII9III9IIIIIIIIIIIIIIIIIII9  
 @ read 108111                                      
 AGCGTGACGGTGTAACGCCCGCCTTTTGATGACTGGGTTTCAAAGAAACG  
 +position: -3336785                                 
 999999999999999IIIIIIII9IIIIIIIII9II99999999999999  
 ```
- `--trim endsfrac (--percent p)`: first trim the ends as in the `ends` option. Accept 
   the trimmed read if the number of low quality nucleotides does not exceed `p%` (default
   `p = 5`).
   Redirect the read to `*_lowq.fq.gz` otherwise.
    Examples (-q 27 [<]): 
 ```
 @ read 1081133
 GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA  GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGG
 +position: -144740                                  +position: -144740 TRIMQ:0:48
 IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 3435                                         @ read  3435 
 CAGTTCTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTGACCCA  CTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTG
 +position: -19377                                   +position: -19377  TRIMQ:5:44
 99999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99999  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 108110                                       
 CAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAAC  
 +position: -4142524                                 
 9IIIIIIIIIIIIIIIIIII9IIII9III9IIIIIIIIIIIIIIIIIII9  
 @ read 108111                                      
 AGCGTGACGGTGTAACGCCCGCCTTTTGATGACTGGGTTTCAAAGAAACG  
 +position: -3336785                                 
 999999999999999IIIIIIII9IIIIIIIII9II99999999999999  
 ```

- `--trim global --global n1 n2`: cut all read globally `n1` nucleotides from 
   the left and `n2` from the right. 

**Note:** qualities are evaluated assuming the reads to follow the
L - Illumina 1.8+ Phred+33, convention, see [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).
Adjust the values for a different convention.            


#### N trimming 

We allow for the following options: 

- `--trimN no`(or flag absent): Nothing is done to the reads containing N's. 
- `--trimN all`: All reads containing at least one N are redirected to `*_NNNN.fq.gz` 
- `--trimN ends`: N's are trimmed if found at the ends, left "as is" otherwise. Example: 
```
@ read 1037                                        @ read 1037
NNTCGACACAAGAAAATGCGCCAATTTTGAGCCAGACCCCAGTTACGCNN TCGACACAAGAAAATGCGCCAATTTTGAGCCAGACCCCAGTTACGC
+position: -2044610                                +position: -2044610  TRIMN:2:47
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 1038                                        
AAAAAACGGTCGTGGTNGGCTACTGTTATTAAAGCGTTGGCTACAAAAAG 
+position: 1361068                                 
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII 
@ read 1039                                        
CCAACGGACNATGGAAGTGCTNCGGCTGGNGTTTTTATCCTCCGGCATTC 
+position: -4282223                                
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII 
@ read 1043                                        @ read 1043 
AATCTGAAACGAAANCTGACAGCGNNCCCCGCTTCTGACAAAATAGGCGN AATCTGAAACGAAANCTGACAGCGNNCCCCGCTTCTGACAAAATAGGCG
+position: -4685876                                +position: -4685876  TRIMN:0:48
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```


- `--trimN strip`: Obtain the largest N free subsequence of the read. Accept it 
   if is longer than the half of the original read length, redirect it to 
   `*_NNNN.fq.gz` otherwise. Example:
```
@ read 1037                                        @ read 1037
NNTCGACACAAGAAAATGCGCCAATTTTGAGCCAGACCCCAGTTACGCNN NNTCGACACAAGAAAATGCGCCAATTTTGAGCCAGACCCCAGTTACGCNN
+position: -2044610                                +position: -2044610  TRIMN:2:47 
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 1038                                        @ read 1038 
AAAAAACGGTCGTGGTNGGCTACTGTTATTAAAGCGTTGGCTACAAAAAG GGCTACTGTTATTAAAGCGTTGGCTACAAAAAG
+position: 1361068                                 +position: 1361068  TRIMN:17:49
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 1039 
CCAACGGACNATGGAAGTGCTNCGGCTGGNGTTTTTATCCTCCGGCATTC 
+position: -4282223                                
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII 
@ read 1043                                         
AATCTGAAACGAAANCTGACAGCGNNCCCCGCTTCTGACAAAATAGGCGN
+position: -4685876 
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```


## Test/examples

  In folder examples, the code is tested in the following way:

1. `./create_fq.sh` generates `EColi_rRNA.fq` that contains:
   * 2e5 reads of length 50 from `EColi_genome.fa` with NO errors.
   * 5e4 reads of length 50 from `rRNA_CRUnit.fa` with NO errrors (rRNA contaminations).
   * Artificially generated reads with low quality score (see script).
   * Artificially generated reads with Ns (see script).
2. The code was tested with flags: 
   `../bin/trimFilter -l 50 --ifa rRNA_CRUnit.fa --ifq EColi_rRNA.fq.gz --tree --trimQ ends --trimN strip`
   i.e., we check for contaminations from rRNA, trim reads with lowQ at the end, and 
   strip reads containing N's. Output is stored in `out*.fq.gz` files. 
3. One can run the example in `./run_example.sh` that will generate the output again 
   with the prefix `tmp`. 
4. With this set up, it is possible to run further customized tests. 


 
## Limitations
  
  The tree option search  limits us to `~5MB` fasta files. 
  Otherwise the tree occupies too much memory. SA search
  is work in progress. 


## Contributors

Paula Pérez Rubio 

## License

GPL v3 (see LICENSE.txt)
