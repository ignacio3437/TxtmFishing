#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import subprocess
from multiprocessing import Pool
import sys
from pathlib import Path

def seq_matcher(source_dir, key_dir, out_path):
    """Will save fasta files in the outpath of 
    a subset of source_dir of the sequences that were
    present in key_dir. 
    IE all of the seqs with the same name as the seqs
     from key_dir will be extracted
    from source_dir and saved to the out_path
    """ 
    loci_path_dict = {}
    key_p = Path(key_dir)
    source_p = Path(source_dir)

    for kfile in [x for x in key_p.iterdir() if x.is_file()]:
        source_file = source_p.joinpath(kfile.name)
        source_dict = SeqIO.to_dict(SeqIO.parse(str(source_file),"fasta"))
        to_write = []
        for record in SeqIO.parse(str(kfile),"fasta"):
            to_write.append(source_dict[record.id])
        with open(str(Path(out_path).joinpath(kfile.name)),'w') as out_handle:
            SeqIO.write(to_write,out_handle,"fasta")
    return





    #     file_prefix = os.path.splitext(ifile)[0]
    #     try:
    #         loci_path_dict[file_prefix].append(os.path.join(idir,ifile))
    #     except KeyError:
    #         loci_path_dict[file_prefix]=[os.path.join(idir,ifile)]
    # for loci in loci_path_dict:
    #     loci_list = loci_path_dict[loci]
    #     for loci_path in loci_list:
    #         loci_out_path = os.path.join(out_path,f"{loci}.fa")
    #         gene_records = []
    #         if isafasta(loci_path):
    #             for record in SeqIO.parse(loci_path, "fasta"):
    #                 gene_records.append(record)
    #             for record in gene_records:
    #                 record.id = record.id.split('-')[0]
    #                 record.description = ''
    #             with open(loci_out_path, 'a') as out_handle:
    #                 SeqIO.write(gene_records, out_handle, "fasta")
    #         else:
    #             print('did not open:'+ loci_path)



source_dir = '/Users/josec/Desktop/Crinoid_capture/MarAA/Results/MarAA_bygene/MarAA_Dna_gen/MarAA_raw_DNA'
key_dir = '/Users/josec/Desktop/Crinoid_capture/MarAA/Results/MarAA_bygene/MarAA_Dna_gen/MarAA_nostop_aln'
out_path = '/Users/josec/Desktop/Crinoid_capture/MarAA/Results/MarAA_bygene/MarAA_Dna_gen/MarAA_DNA'

seq_matcher(source_dir,key_dir,out_path)