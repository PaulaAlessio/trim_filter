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


```
Usage: trimFilter --ifq INPUT.fq -l READLEN [-q MINQ | --ifa INPUT.fa | --idx IDX_FILE | 
                  -o O_PREFIX | (--tree|-t) | --trimQ | --trimN [all|ends] | (--help|-h) ]
Reads in a fq file (gz, bz2, z formats also accepted) and removes: 
  * Low quality reads (trim reads passing the --trim option),
  * Reads containing N base callings,
Outputs 4 [O_PREFIX]_fq.gz files containing: "good" reads, discarded low Q reads,
discarded reads containing N's, discarded contaminations.
Options:
 -h, --help   prints help dialog.
 -l           Read length: length of the reads. Required option.
 -f, --ifq    Fastq input file [*fq|*fq.gz|*fq.bz2]: Required option.
 -a, --ifa    Fasta input file [*fa|*fa.gz|*fa.bz2]: Optional.
 -x, --idx    Index input file: Optional.
 -k, --lmer   Lmer length: Length of the Lmers to be searched on the contaminations *fa.
              Default: length of the reads.
 -t, --tree   Construct a tree for the search process. If missing, SA constructed.
 -Q, --trimQ  trims lowQ reads if present. Discards all reads with lowQ bases if absent.
 -N, --trimN  (all, ends) trims reads containing N's if all, only trims at the ends 
              if ends and discards all reads with N's if absent.
 -q           Minimum quality allowed. Optional option. Default 27 .
 -o           Output prefix (containing the path). Optional option. Default ./out.
```

## Output description

- `O_PREFIX_good.fq.gz`: contains reads that passed all filters.
- `O_PREFIX_NNNN.fq.gz`: contains reads with *too many* N's. 
- `O_PREFIX_lowQ.fq.gz`: contains reads discarded due to low quality issues.
- `O_PREFIX_cont.fq.gz`: contains contamination reads (matching the fasta sequence). 


## Filters

   Filters were applied in the following fashion: 

####**LowQ**

 **Modify the CODE: 3 options: no all value of lowQ**

  Look for lowQ (below minQ) base callings at the beginning and at the end of the read. 
  Trim them at both ends until the quality is above the threshold. Keep the read 
  and annotate in the fourth line where the read has been trimmed (starting to count
  from 0) if:
    - the read has no further qualities below the threshold. 
    - The length of the remaining part is larger than the half of the original read length. 

  Examples (-q 27 [<]): 
```
@ read 3435                                          @ read  3435 
CAGTTCTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTGACCCA   CTTTGGTGTAGATGGAGCAGGACGGGATACCCATACCGTG
+position: -19377                                    +position: -19377  TRIMQ:5:44
99999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99999   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 108110                                        @ read 108110 
CAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAAC   AGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAA
+position: -4142524                                  +position: -4142524  TRIMQ:1:48
9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 108111                                       
AGCGTGACGGTGTAACGCCCGCCTTTTGATGACTGGGTTTCAAAGAAACG   
+position: -3336785                                  
99IIIIIIIIIIII9IIIIIIII9IIIIIIIII9IIIIIIIIIIIIII99   
@ read 108112                                       
ACCGGTAATGATACCCGCCAGTACACCCATGTTTATTTCTGGGTTGATGG   GTAATGATACCCGCCAGTACACCCATGTTTATTTCTGGGTTGATGG
+position: -4548423                                  +position: -4548423 TRIMQ:4:49
9999IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 1081133                                       @ read 1081133 
GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGGA   GTACATTGATGAGCAACATGGGGCTTGAACTGGCGCTGAAACAGTTAGG
+position: -144740                                   +position: -144740 TRIMQ:0:48
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@ read 108114                                       
CTTAACCTGAACATTGAAGAGTTTAAAGGCGAGCAGTTACGCGTGACGGG   
+position: 5075572                                      
IIIIIIIIIIIIIIIIIIIII9999999IIIIIIIIIIIIIIIIIIIIII     
@ read 108116                                       
TCAGGGTGCAGCGGTAGCGCCAGTTTCCAGTCCATCGCCATTTTGTAGAC   
+position: 1750468                                 
9IIIIIIIIIII9IIIIIIIIIII9IIIIIIII9IIIIIIIIIIIIIII9 
```


**Note:** qualities are evaluated assuming the reads to follow the
L - Illumina 1.8+ Phred+33, convention, see [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).
Adjust the values for a different convention.            

#### **N trimming  **

We allow for the following options: 

- `no`(default): If the flag is not 
- `ends`
- `all`
- `strip`

**MODIFY THE CODE ACCORDINGLY**

## Test/examples

  In folder examples, the code is tested in the following way: 

**FINISH eRITING**

 
## Limitations
  
  Only the tree option search is implemented. This limits us 
to `~5MB` fasta files. Otherwise the tree occupies too much memory.


## Contributors

Paula Perez Rubio 

## License

GPL v3 (see LICENSE.txt)
