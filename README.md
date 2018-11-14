# TxtmFishing

TxtmFishing pulls target protein coding genes from transcriptomes and removes portions of the gene that were not targeted. An optimal blast stringency is chosen by iterating through blast e-values until a minimum number of transcriptomes with hits is reached. The minimum number of transcriptomes and a list of e-values to test are specified in the parameter file.

Output files are formated to be used as a reference 'target file' for HybPiper.

Usage:
```
python3 TxtmFishing.py --param parameters.txt 
```



## Input Files

1) Transcriptome files formated as a genbank query databases in a folder.
2) Full CDS of target genes as AA.fa
3) Concatenated target exons of genes as AA.fa
4) Text file of transcriptome names.
5) Text file of a list of the gene names.
6) Parameter file. See Param-default.txt

## Dependencies
+ tBlastx
+ Exonerate
+ Mafft

Python3 Libraries:
+ Pool
+ BioPython
+ regex
