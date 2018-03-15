#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import subprocess
from multiprocessing import Pool
import sys
from os.path import basename


def dir_CATer(dir1, dir2, out_path):
    """Combines the files in two directories.
    If the same fasta is in both directories, they are concatenated.
    """ 
    loci_path_dict = {}
    directories=[dir1,dir2]
    for idir in directories:
        for ifile in os.listdir(idir):
            file_prefix = os.path.splitext(ifile)[0]
            try:
                loci_path_dict[file_prefix].append(os.path.join(idir,ifile))
            except KeyError:
                loci_path_dict[file_prefix]=[os.path.join(idir,ifile)]
    for loci in loci_path_dict:
        loci_list = loci_path_dict[loci]
        for loci_path in loci_list:
            loci_out_path = os.path.join(out_path,f"{loci}.fa")
            gene_records = []
            if isafasta(loci_path):
                for record in SeqIO.parse(loci_path, "fasta"):
                    gene_records.append(record)
                for record in gene_records:
                    record.id = record.id.split('-')[0]
                    record.description = ''
                with open(loci_out_path, 'a') as out_handle:
                    SeqIO.write(gene_records, out_handle, "fasta")
            else:
                print('did not open:'+ loci_path)