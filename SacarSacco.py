#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet


"""
Removes all seqs with rem_name from all of the fasta files in in_dir. 
Usage: Sacalo.py --in ind_dir --name rem_name --out out_dir 
"""

in_dir = "/Users/josec/Desktop/Crinoid_capture/MarAA/Results/alignments/MarAA_DNA_aln/nodem"
rem_name = "PSATHY"
out_dir = "/Users/josec/Desktop/Crinoid_capture/MarAA/Results/alignments/MarAA_DNA_aln/trim"

try:
    os.mkdir(out_dir)
except:
    pass

for in_file_path in os.listdir(in_dir):
    if in_file_path.lower().endswith(('.fa', '.fasta', '.fna')):
        seq_records = SeqIO.parse(os.path.join(in_dir,in_file_path), "fasta")
        with open(f'{out_dir}/{in_file_path}', 'w') as outhandle:
            for record in seq_records:
                if rem_name in record.id:
                    pass
                else:
                    SeqIO.write(record, outhandle, "fasta")