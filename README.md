# TxtmFishing
Pulls orthologous genes from transcriptomes and trims them to target exons. Output files are formated to be used as a reference for HybPiper.

Input Files:
1) Transcriptome files formated as a genbank query database.
2) Full target genes as aa.fa
3) Targeted exons of genes as dna.fa
4) Text file of list of transcriptomes
5) Text file of a list of the gene names
6) Parameter file. See example

Required:

Python3, Pool, BioPython, regex

tBlastx, Exonerate, Mafft
