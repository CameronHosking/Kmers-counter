# Welcome to the Kmers-counter wiki!
# Arguments
Arguments are supplied as either flags or key value 
pairs separated with an equals sign.

## Required arguments

**-input=<filename>** 
This specifies the input file which must be in fasta format or just a string of characters

**-K=<integer>** 
This specifies the length of kmers to count, currently supports values between 1 and 15 inclusive


## Optional arguments

**-outputZeros**  
If this flag is present than kmers with counts of 
zero are output otherwise they are skipped

**-split** 
If this flag is present then kmer counts between id tags (which are defined as lines starting with '>') 
are counted separately and output into separate files. Each file is headed by the id tag. the name of each file is simply output+position of tag within file

## Output
The output is one or more files containing the frequency of kmers occurring in the input file.
The output is a tab seperated file in which the Kmers are listed alphabetically in column one and their associated frequencies are in column two.

## Examples
    ./KmersCounter -input=human_g1k_v37.fasta -K=8  
This will read in the file `human_g1k_v37.fasta` and 
output the frequency of all the different sequences 
consisting of [ACGT] of length 8 to one file.

##

    ./KmersCounter -input=human_g1k_v37.fasta -K=8 -outputZeros
This will do the same as above but will also output sequences that have a frequency of zero.

##

    ./KmersCounter -input=test.txt -K=3
where test.txt is just `CAGTAGCNGATCAGGTAGCCGATNAGCANTAGGANNNNNNNNNNN`
produces the output file output54.txt.  

    AGC	0.130435
    AGG	0.0869565
    AGT	0.0434783
    ATC	0.0434783
    CAG	0.0869565  
    CCG	0.0434783  
    CGA	0.0434783  
    GAT	0.0869565  
    GCA	0.0434783  
    GCC	0.0434783  
    GGA	0.0434783  
    GGT	0.0434783  
    GTA	0.0869565  
    TAG	0.130435  
    TCA	0.0434783
