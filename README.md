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
                  -o O_PREFIX | (--tree|-t) | --trimQ | --trimN [all|ends] | (--help|-h) ]
Reads in a fq file (gz, bz2, z formats also accepted) and removes: 
  * Low quality reads (trim reads passing the --trim option),
  * Reads containing N base callings,
Outputs 4 [O_PREFIX]_fq.gz files containing: "good" reads, discarded low Q reads,
discarded reads containing N's, discarded contaminations.
 -h, --help   prints help dialog.
 -l           Read length: length of the reads. Required option.
 -f, --ifq    Fastq input file [*fq|*fq.gz|*fq.bz2]: Required option.
 -a, --ifa    Fasta input file [*fa|*fa.gz|*fa.bz2]: Optional.
 -x, --idx    Index input file: Optional.
 -k, --lmer   Lmer length: Length of the Lmers to be searched on the contaminations *fa.
              Default: length of the reads.
 -t, --tree   Construct a tree for the search process. If missing, SA constructed.
 -Q, --trimQ  no (or no flag): does nothing to low quality reads,
              ends: trims the ends of the read if their quality is below the threshold -q,
              all:  removes all reads containing at least one low quality nucleotide.
 -N, --trimN  no (or no flag): does nothing to reads containing N's,
              all: removes all reads containing N's,
              ends: trims ends of reads with N's,
              strips: looks for the largest substring containing no N's.
 -q           Minimum quality allowed. Optional option. Default 27 .
 -o           Output prefix (containing the path). Optional option. Default ./out.
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

#### LowQ

- `--trimQ no` or flag absent: nothing is done to the reads with low quality.
- `--trimQ all`: all reads containing at least one low quality nucleotide are
  redirected to  `*_lowq.fq.gz`
- `--trimQ ends`: look for low quality (below minQ) base callings at the beginning and at the end of the read. 
  Trim them at both ends until the quality is above the threshold. Keep the read in `*_good.fq.gz`
  and annotate in the fourth line where the read has been trimmed (starting to count
  from 0) if:
    - the read has no further qualities below the threshold. 
    - The length of the remaining part is larger than the half of the original read length. 

    Examples (-q 27 [<]): 
 ```
 @ read 3435                                         @ read  3435 
 CAGTTCTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTGACCCA  CTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTG
 +position: -19377                                   +position: -19377  TRIMQ:5:44
 99999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99999  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 108110                                       @ read 108110 
 CAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAAC  AGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAA
 +position: -4142524                                 +position: -4142524  TRIMQ:1:48
 9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 108111                                      
 AGCGTGACGGTGTAACGCCCGCCTTTTGATGACTGGGTTTCAAAGAAACG  
 +position: -3336785                                 
 99IIIIIIIIIIII9IIIIIIII9IIIIIIIII9IIIIIIIIIIIIII99  
 @ read 108112                                      
 ACCGGTAATGATACCCGCCAGTACACCCATGTTTATTTCTGGGTTGATGG  GTAATGATACCCGCCAGTACACCCATGTTTATTTCTGGGTTGATGG
 +position: -4548423                                 +position: -4548423 TRIMQ:4:49
 9999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 1081133                                      @ read 1081133 
 GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA  GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGG
 +position: -144740                                  +position: -144740 TRIMQ:0:48
 IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 @ read 108114                                      
 CTTAACCTGAACATTGAAGAGTTTAAAGGCGAGCAGTTACGCGTGACGGG  
 +position: 5075572                                     
 IIIIIIIIIIIIIIIIIIIII9999999IIIIIIIIIIIIIIIIIIIIII    
 @ read 108116                                      
 TCAGGGTGCAGCGGTAGCGCCAGTTTCCAGTCCATCGCCATTTTGTAGAC  
 +position: 1750468                                
 9IIIIIIIIIII9IIIIIIIIIII9IIIIIIII9IIIIIIIIIIIIIII9 
 ```
 Redirect the read to `*_lowq.fq.gz` otherwise.

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
+position: -4685876                                +position: -4685876  TRIMN:0-48
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
  
  Only the tree option search is implemented. This limits us 
to `~5MB` fasta files. Otherwise the tree occupies too much memory.


## Contributors

Paula Pérez Rubio 

## License

GPL v3 (see LICENSE.txt)
