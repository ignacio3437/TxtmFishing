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
from pathlib import Path


def filelines_to_list(file):
    """Makes a list of each new line in a file. """
    with open(file, 'rU') as file_handle:
        file_list = [ind.rstrip() for ind in file_handle.readlines()]
    return file_list

def run_command(command,command_out_path):
    command_output = subprocess.getstatusoutput(command)
    with open(command_out_path, "w") as out_handle:
        out_handle.write(command_output[1])
        if command_output[0]:
            print(f"{command} Failed")
            sys.exit(1)
    return

def out_dir_maker(project_path):
    """Creates output directories for the stages of the pipeline. """
    dirs_to_make = [f"{project_path}/Blast_results", f"{project_path}/Txtms_TrimtoTargets_alignments",
                    f"{project_path}/Blast_by_gene", f"{project_path}/Exonerate_out", f"{project_path}/Exonerate_clean"]
    for directory in dirs_to_make:
        try:
            os.mkdir(directory)
        except IOError:
            pass
    return dirs_to_make


def blast_run(evalue, threads, in_path, query_file, out_path, org_file):
    """Run tblastn on each file in in_path with set evalue cutoff on x threads.
    Only saves top blast hit. Infile must end with '.nt'.
    Consecutive runs will overwrite the blast_output file.
    """
    filenames = os.listdir(in_path)
    orglist = filelines_to_list(org_file)
    for file in filenames:
        if file in orglist:
            filepath = str(Path(in_path)/file)
            # file_prefix = file.replace('.Trinity.annotated.nt', '')
            command = ' '.join((f"tblastn", 
            f"-query {query_file} -db {filepath}",
            f"-evalue {evalue} -num_threads {threads}",
            f"-outfmt '6 sseqid qseqid' -num_alignments 1"))
            command_out_path = f"{out_path}/{file}.blast.txt"
            run_command(command, command_out_path)
    return


def seq_fetcher(list_dir, list_file_endings, seq_dir, org_file):
    """Loads the blast results of each organism and writes the results by gene.
    Appends to existing files if present."""
    orglist = filelines_to_list(org_file)
    for org in orglist:
        if len(org)>3: #Make sure the line is not blank
            hit_records = []
            record_dict = SeqIO.to_dict(SeqIO.parse(
                f"{seq_dir}/{org}", "fasta"))
            tofetcha = filelines_to_list(f"{list_dir}/{org}{list_file_endings}")
            tofetch = [tf for tf in tofetcha if 'Warning' not in tf]
            for fetch_id in tofetch:
                record_id, gene_id = fetch_id.split()
                record = record_dict[record_id]
                record_i = SeqRecord(
                    Seq(f"{record.seq}", SingleLetterAlphabet()), id = gene_id)
                hit_records.append(record_i)
            out_file_path=f'{list_dir}/{org}.fa'
            with open(out_file_path, 'w') as outhandle:
                SeqIO.write(hit_records, outhandle, "fasta")
            fasta_deduper(fasta_file=out_file_path)

    return


def cat_by_gene(org_file, loci_list, in_path, file_ending, out_path, cutoff):
    """This concatenates sequence files for all organisms by loci.
    Sequence files must have fasta heading of '>org-loci'
    The name of the alignments with fewer sequences than 'cutoff' are returned in check_list."""
    check_list = []
    org_list = filelines_to_list(org_file)
    org_dict = {}
    for org in org_list:
        if len(org)>3: #Make sure the line is not blank
            exonerate_dict = SeqIO.to_dict(SeqIO.parse(
                f"{in_path}/{org}{file_ending}", "fasta"))
            org_dict[org] = exonerate_dict
    spurgenes = loci_list
    # print(f"Number of genes: {len(spurgenes)}")
    # print(f"Number of transcriptomes: {len(org_list)}")
    for spurgene in spurgenes:
        out_file_path = f"{out_path}/{spurgene}.fa"
        with open(out_file_path, 'w') as out_handle:
            for org in org_list:
                record_id = f"{org}-{spurgene}"
                try:
                    record_i = SeqRecord(Seq(f"{org_dict[org][spurgene].seq}", SingleLetterAlphabet()), id=record_id, description="")
                    SeqIO.write(record_i, out_handle, "fasta")
                except KeyError:
                    pass
        with open(out_file_path, 'r') as out_handle:
            records = list(SeqIO.parse(out_handle, "fasta"))
            if len(records) < cutoff:
                check_list.append(spurgene)
    return check_list


def exoneratetor_bygene(query_file, in_path, in_file_ending, out_path, loci_file, num_threads):
    """Makes a list of commands to send to exonerator() and parallizes the run. """
    loclist = filelines_to_list(loci_file)
    command_args = []
    for loc in loclist:
        in_file = f"{in_path}/{loc}{in_file_ending}"
        command_args.append(
            [query_file, in_file, in_file_ending, out_path, loc])
    p = Pool(num_threads)
    p.map(exonerator, command_args)
    return


def exonerator(command_args):
    """Run exonerate. This is a function to allow parallelization with multiprocessing.Pool()"""
    query_file, in_file, in_file_ending, out_path, loc = command_args
    command_out_path=f"{out_path}/{loc}{in_file_ending}"
    command= ' '.join((
        f"exonerate --model protein2genome",
        f"-q {query_file} -t {in_file}",
        f"-Q protein -T dna",
        f"--showvulgar F --showalignment F --verbose 0 --fsmmemory 20G --ryo '>%ti %qi\\n%tcs\\n'"))
    # print(command)
    # print(command_out_path)
    run_command(command,command_out_path)
    return


def exonerateout_cleaner(exonerateout_path, file_ending, loci_file, exonerateclean_path):
    """This parses the exonerate output from exonerator so that the query hit and the target hit are the same loci."""
    loclist = filelines_to_list(loci_file)
    fasta_paths_to_dedupe = []
    for loc in loclist:
        to_add = []
        in_file = f"{exonerateout_path}/{loc}{file_ending}"
        records = SeqIO.parse(in_file, "fasta")
        for record in records:
            qloci = record.description.split()[1]
            tloci = record.id.split('-')[1]
            if qloci == tloci:
                to_add.append(record)
        out_path = f"{exonerateclean_path}/{loc}{file_ending}"
        to_write = []
        for record in to_add:
            org = record.id
            record_name = f"{org.split('.')[0]}"
            record.id = record_name
            to_write.append(record)
        with open(out_path, 'w') as out_handle:
            SeqIO.write(to_write, out_handle, 'fasta')
            fasta_paths_to_dedupe.append(os.path.abspath(out_path))
    return fasta_paths_to_dedupe


def fasta_deduper(fasta_file):
    """Rewrites over fasta file with only one seq per name. Keeps the longest sequence. Removes seq description."""
    u_records = []
    u_record_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        record.description = ''
        if record.id in u_record_dict:
            if u_record_dict[record.id][0] < len(record.seq):
                u_record_dict[record.id] = [len(record.seq), record]
            else:
                pass
        else:
            u_record_dict[record.id] = [len(record.seq), record]
    u_records = [u_record_dict[record][1] for record,val in u_record_dict.items()]
    with open(fasta_file, 'w') as out_handle:
        SeqIO.write(u_records, out_handle, "fasta")
    return len(u_record_dict)


def fasta_subseter(subset_list, in_fasta_path, out_fasta_path):
    """Given a list of fasta headers, it will extract the seqs in that list from in_fasta_path
    and write them to out_fasta_path
    """
    input_dict = SeqIO.to_dict(SeqIO.parse(in_fasta_path, "fasta"))
    out_records = []
    for name in subset_list:
        out_records.append(input_dict[name])
    with open(out_fasta_path, "w") as out_handle:
        SeqIO.write(out_records, out_handle, 'fasta')
    return


def maffter(fasta_file_path, out_file_path, num_threads):
    """ Runs and auto mafft alignemnt on a fasta file.
    Settings changed to produce a less gappy alignment because we are dealing with coding data. """
    command = f" mafft --oldgenafpair --leavegappyregion --thread {num_threads} --op 2 --quiet {fasta_file_path}"
    command_out_path = out_file_path
    run_command(command,command_out_path)
    return

def isafasta(fasta_to_check_path):
    """Returns True if file is parsable as a fasta by BioPython"""
    records = SeqIO.parse(fasta_to_check_path,"fasta") 
    try:
        return any(records)
    except:
        return False

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
    return


def param_reader(paramfile_path):
    """Reads a parameter file with the following:
    Project project_path
    Number of threads
    Minimum number of sequences for a given e_vals blast cutoff
    The blast evalue cutoffs to test in the order to test. Formated as python list.
    Filename for blast query file
    Filename for the list of transcriptomes
    Filename for the list of loci
    Filename of the exonerate query
    Path of Transcriptomes
    """
    plist = filelines_to_list(paramfile_path)
    p2_list = [p.split('#')[0].strip() for p in plist]
    param_list = []
    param_list.append(p2_list[0])
    for x in p2_list[1:3]:
        param_list.append(int(x))
    e_vals = [x.replace(' ','') for x in p2_list[3].split(',')]
    param_list.append(e_vals)
    for x in p2_list[4:]:
        param_list.append(os.path.join(p2_list[0],x)) 
    return param_list

def txtm_fishing_pipe(param_list):
    """Run pipeline in Python.
    This will iterate the blast searches over the e_vals.
    If a gene has less than min_num_seqs for a given evalue, the evalue is lowered and the search is repeated.
    The results from all of the blast searches are writen to a file and fed into exonerate.
    Exonerate removes exons that we did not target (ie any exons not present in the exonerate_query seq) and regions that do not align with the reference. 
    """
    project_path, num_threads, min_num_seqs, e_vals, blast_query_path, org_list_path, loci_list_path, exonerate_query_path, txtm_folder = param_list
    dirs_to_make = out_dir_maker(project_path)
    blast_results, mafft_alignments, blast_by_gene, exonerate_out, exonerate_clean = dirs_to_make
    to_blast_list = filelines_to_list(loci_list_path)
    # Start the pipeline:
    for e_val in e_vals:
        blast_run(evalue=e_val, threads=num_threads, in_path=txtm_folder, query_file=blast_query_path, out_path=blast_results, org_file=org_list_path)
        seq_fetcher(list_dir=blast_results, list_file_endings='.blast.txt', seq_dir=txtm_folder, org_file=org_list_path)
        check_list = cat_by_gene(org_list_path, to_blast_list, blast_results, '.fa', blast_by_gene, min_num_seqs)
        #Print the status of the blast searches
        if len(check_list) < 20:
            print(
                f"The following loci had no hits for {e_val}:{' '.join(check_list)}")
        elif not check_list:
            print("All loci had hits at {e_val}")
        else:
            print(
                f"There were {len(check_list)} genes with < {min_num_seqs} hits for {e_val}")
        if not check_list:
            break
        else:
            to_blast_list = check_list
            new_query_path = f'{project_path}/SpurGenes_aa_{e_val}.fasta'
            fasta_subseter(to_blast_list, blast_query_path, new_query_path)
            blast_query_path = new_query_path
    exoneratetor_bygene(exonerate_query_path, blast_by_gene, '.fa', exonerate_out, loci_list_path, num_threads)
    fasta_paths_to_dedupe = exonerateout_cleaner(exonerate_out, '.fa', loci_list_path, exonerate_clean)
    for fasta_file in fasta_paths_to_dedupe:
        num_seqs = fasta_deduper(fasta_file)
        if num_seqs > 1:
            maffter(fasta_file, os.path.join(mafft_alignments, os.path.basename(fasta_file)), num_threads)
        else:
            pass
    return exonerate_clean

def wrap_up(directory_to_wrap_up,outfile_path=None):
    """Wraps up a directory by collating all of the fasta sequences in that directory into a single file.
    outfile_path is the baitfile to be used in HybPiper.
    """
    out_records = []
    for ifile in os.listdir(directory_to_wrap_up):
        gene = ifile.split('.')[0]
        for irecord in SeqIO.parse(os.path.join(directory_to_wrap_up,ifile),"fasta"):
            taxa = irecord.id
            irecord.id = f"{taxa}-{gene}"
            irecord.description = ''
            out_records.append(irecord)
    num_records = (len(out_records))
    SeqIO.write(out_records, outfile_path, "fasta")
    return num_records

def main():
    # Run the pipeline
    args = sys.argv[1:]
    # args = ["--param", "/Users/josec/Desktop/git_repos/TxtmFishing/Mar-param.txt"] # For Testing
    usage = 'usage: TxtmFishing.py --param parameters.txt'
    if not args:
        print(usage)
        sys.exit(1)
    if args[0] == '--param':
        paramfile_path = args[1]
    else:
        print(usage)
        sys.exit(1)
    param_list = param_reader(paramfile_path)
    dir2wrap = txtm_fishing_pipe(param_list)
    baitfile_path = f"{param_list[0]}/baitfile.fasta"
    num_records = wrap_up(dir2wrap,outfile_path=baitfile_path)
    print(f"{num_records} seqs written to baitfile:{baitfile_path}")
    print('\a')
if __name__ == "__main__":
    main()
