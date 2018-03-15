#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet


"""
Removes all seqs with stop codons from all of the fasta files in in_dir. 
"""

in_dir = "/Users/josec/Desktop/Crinoid_capture/MarAA/Results/MarAA_bygene/aa/"
stop_codon = "*"
out_dir = "/Users/josec/Desktop/Crinoid_capture/MarAA/Results/MarAA_bygene/aa_nostop/"

try:
    os.mkdir(out_dir)
except:
    pass

for in_file_path in os.listdir(in_dir):
    if in_file_path.lower().endswith(('.fa', '.fasta', '.fna')):
        seq_records = SeqIO.parse(os.path.join(in_dir,in_file_path), "fasta")
        with open(f'{out_dir}/{in_file_path}', 'w') as outhandle:
            for record in seq_records:
                if stop_codon in str(record.seq):
                    pass
                else:
                    SeqIO.write(record, outhandle, "fasta")