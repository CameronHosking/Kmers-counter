# Welcome to the Kmers-counter wiki!
# Arguments
Arguments are supplied as either flags or key value 
pairs separated with an equals sign.

## Required arguments

**-K=<integer>** 
This specifies the length of kmers to count, currently supports values between 1 and 15 inclusive


## Optional arguments

**-outputZeros**  
If this flag is present than kmers with counts of 
zero are output otherwise they are skipped

**-split** 
If this flag is present then kmer counts between id tags (which are defined as lines starting with '>') 
are counted separately and output into separate files. Each file is headed by the id tag.

**-input=<filename>** 
This specifies the input file which must be in fasta format or just a string of characters. If this option is not specified then the sequence will be read from stdin.

**-output=<filename>** 
This specifies the output file(s). If the split flag is present then the files will be numbered from 1 to the number of id lines present. to select where the number appears in the filename add a %, the first occurrence of which will be replaced by the file number. If this option is not specified output will be to stdout.

## Output
The output is tab seperated in which the Kmers are listed alphabetically in column one and their associated frequencies are in column two.

## Examples
    ./KmersCounter -K4

This will read a sequence from stdin and output the kmer frequencies to stdout only showing thos kmers that appeared at least once

##

    ./KmersCounter -input=human_g1k_v37.fasta -K=8 -output="out.txt"  
This will read in the file `human_g1k_v37.fasta` and 
output the frequency of all the different sequences 
consisting of [ACGT] of length 8 to one file "out.txt".

##

    /KmersCounter -K=8 -output=out -split 
This will from stdin and output the frequency of all the different sequences consisting of [ACGT] of length 8 seperated by id to 1 file per id: "out1", "out2", "out3", etc.

##

    /KmersCounter -K=8 -output="out%.txt" -split 
This will be the same as above except the files will be named "out1.txt", "out2.txt", "out3.txt", etc.

##

    ./KmersCounter -input=human_g1k_v37.fasta -K=8 -outputZeros
This will will also output sequences that have a frequency of zero.

##

    ./KmersCounter -input=test.txt -K=3
where test.txt is just `CAGTAGCNGATCAGGTAGCCGATNAGCANTAGGANNNNNNNNNNN`
produces the output.  

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