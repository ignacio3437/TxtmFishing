#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import subprocess
from multiprocessing import Pool
import sys


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
    dirs_to_make = [f"{project_path}/Blast_databases", f"{project_path}/Blast_results", f"{project_path}/Txtms_TrimtoTargets_alignments",
                    f"{project_path}/Blast_by_gene", f"{project_path}/Exonerate_out", f"{project_path}/Exonerate_clean"]
    for directory in dirs_to_make:
        try:
            os.mkdir(directory)
        except IOError:
            pass
    return dirs_to_make


def blast_run(evalue, threads, in_path, query_file, out_path):
    """Run tblastn on each file in in_path with set evalue cutoff on x threads.
    Only saves top blast hit. Infile must end with '.nt'.
    Consecutive runs will overwrite the blast_output file.
    """
    filenames = os.listdir(in_path)
    for file in filenames:
        if file.endswith(".nt"):
            file_prefix = file.replace('.Trinity.annotated.nt', '')
            command = ' '.join((f"tblastn", 
            f"-query {query_file} -db {in_path}/{file}",
            f"-evalue {evalue} -num_threads {threads}",
            f"-outfmt '6 sseqid qseqid' -max_target_seqs 1 -max_hsps 1"))
            command_out_path = f"{out_path}/{file_prefix}.blast.txt"
            run_command(command, command_out_path)
    return


def seq_fetcher(list_dir, list_file_endings, seq_dir, seq_file_ending, org_file):
    """Loads the blast results of each organism and writes the results by gene.
    Appends to existing files if present."""
    orglist = filelines_to_list(org_file)
    for org in orglist:
        hit_records = []
        record_dict = SeqIO.to_dict(SeqIO.parse(
            f"{seq_dir}/{org}{seq_file_ending}", "fasta"))
        tofetch = filelines_to_list(f"{list_dir}/{org}{list_file_endings}")
        for fetch_id in tofetch:
            record_id, gene_id = fetch_id.split()
            record = record_dict[record_id]
            record_i = SeqRecord(
                Seq(f"{record.seq}", SingleLetterAlphabet()), id = gene_id)
            hit_records.append(record_i)
        with open(f'{list_dir}/{org}.fa', 'w') as outhandle:
            SeqIO.write(hit_records, outhandle, "fasta")
    return


def cat_by_gene(org_file, loci_list, in_path, file_ending, out_path, cutoff):
    """This concatenates sequence files for all organisms by loci.
    Sequence files must have fasta heading of '>org-loci'
    The name of the alignments with fewer sequences than 'cutoff' are returned in check_list."""
    check_list = []
    org_list = filelines_to_list(org_file)
    org_dict = {}
    for org in org_list:
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
                    record_i = SeqRecord(Seq(f"{org_dict[org][spurgene].seq}", SingleLetterAlphabet(
                    )), id=record_id, description="")
                    SeqIO.write(record_i, out_handle, "fasta")
                except:
                    pass
        with open(out_file_path, 'r') as out_handle:
            records = list(SeqIO.parse(out_handle, "fasta"))
            if len(records) < cutoff:
                check_list.append(spurgene)
    return check_list


def exoneratetor_bygene(query_file, in_path, in_file_ending, out_path, loci_file, org_file, num_threads):
    """Makes a list of commands to send to exonerator() and parallizes the run. """
    orglist = filelines_to_list(org_file)
    loclist = filelines_to_list(loci_file)
    query_dict = SeqIO.to_dict(SeqIO.parse(query_file, "fasta"))
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
        f"-Q protein -T dna"
        f"--showvulgar F --showalignment F --verbose 0 --fsmmemory 20G --ryo '>%ti %qi\n%tas\n'"))
    run_command(command,command_out_path)
    return


def exonerateout_cleaner(exonerateout_path, file_ending, loci_file, exonerateclean_path):
    """This parses the exonerate output from exonerator so that the query hit and the target hit are the same loci."""
    loclist = filelines_to_list(loci_file)
    fasta_paths_to_dedupe = []
    for loc in loclist:
        to_write = []
        in_file = f"{exonerateout_path}/{loc}{file_ending}"
        records = SeqIO.parse(in_file, "fasta")
        for record in records:
            qloci = record.description.split()[1]
            tloci = record.id.split('-')[1]
            if qloci == tloci:
                to_write.append(record)
        out_path = f"{exonerateclean_path}/{loc}{file_ending}"
        with open(out_path, 'w') as out_handle:
            SeqIO.write(to_write, out_handle, 'fasta')
            fasta_paths_to_dedupe.append(os.path.abspath(out_path))
    return fasta_paths_to_dedupe


def fasta_deduper(fasta_file):
    """Rewrites over fasta file with only one seq per name. Keeps the longest sequence. Removes seq description"""
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
    u_records = [u_record_dict[record][1] for record in u_record_dict.keys()]
    with open(fasta_file, 'w') as out_handle:
        SeqIO.write(u_records, out_handle, "fasta")
    return


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
    """ Runs and auto mafft alignemnt on a fasta file"""
    command = f"mafft --thread {num_threads} --quiet --auto {fasta_file_path}"
    command_out_path = out_file_path
    run_command(command,command_out_path)
    return


def txtm_fishing_pipe(param_list):
    """Run pipeline in Python! Shell scripts are #OldSckool.
    This will iterate the blast searches over the e_vals.
    If a gene has no hits for a given evalue, it is lowered and searched again.
    The results from all of the blast searches are concatenated and fed into exonerate.
    Exonerate picks out the exons that we targeted. 
    """
    project_path, num_threads, min_num_seqs, e_vals, blast_query_path, org_list_path, loci_list_path, exonerate_query_path = param_list
    dirs_to_make = out_dir_maker(project_path)
    blast_databases, blast_results, mafft_alignments, blast_by_gene, exonerate_out, exonerate_clean = dirs_to_make
    to_blast_list = filelines_to_list(loci_list_path)
    # Start the pipeline:
    for e_val in e_vals:
        blast_run(e_val, num_threads, blast_databases,
                  blast_query_path, blast_results)
        seq_fetcher(blast_results, '.blast.txt', blast_databases, '.Trinity.annotated.nt', org_list_path)
        check_list = cat_by_gene(
            org_list_path, to_blast_list, blast_results, '.fa', blast_by_gene, min_num_seqs)
        if len(check_list) < 20:
            print(
                f"The folowing loci had no hits for {e_val}:{' '.join(check_list)}")
        elif len(check_list) == 0:
            print("All loci had hits at {e_val}")
        else:
            print(
                f"There were {len(check_list)} genes with no hits for {e_val}")
        to_blast_list = check_list
        new_query_path = f'{project_path}/SpurGenes_aa_{e_val}.fasta'
        fasta_subseter(to_blast_list, blast_query_path, new_query_path)
        blast_query_path = new_query_path
    exoneratetor_bygene(exonerate_query_path, blast_by_gene, '.fa',
                        exonerate_out, loci_list_path, org_list_path, num_threads)
    fasta_paths_to_dedupe = exonerateout_cleaner(
        exonerate_out, '.fa', loci_list_path, exonerate_clean)
    for fasta_file in fasta_paths_to_dedupe:
        fasta_deduper(fasta_file)
        maffter(fasta_file, os.path.join(mafft_alignments, os.path.basename(fasta_file)), num_threads)
    txtm_tarlen(loci_list_path, org_list_path, exonerate_clean, '/Users/josec/Desktop/Crinoid_capture/Feb5_hybTxCrinoid/pre_Txm_tarlens.txt')
    return


def dir_CATer(dir1, dir2, out_path):
    """Combines the files in two directories.
    If the same fasta is in both directories, they are concatenated.
    """
    set_files = list(set(os.listdir(dir1)+os.listdir(dir2)))
    for file in set_files:
        gene_records = []
        [gene_records.append(record)
         for record in SeqIO.parse(f"{dir1}/{file}", "fasta")]
        [gene_records.append(record)
         for record in SeqIO.parse(f"{dir2}/{file}", "fasta")]
        for record in gene_records:
            record.id = record.id.split('-')[0]
            record.description = ''
        with open(f"{out_path}/{file}", 'w') as out_handle:
            SeqIO.write(gene_records, out_handle, "fasta")
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
    """
    plist = filelines_to_list(paramfile_path)
    p2_list = [p.split('#')[0].strip() for p in plist]
    param_list = []
    param_list.append(p2_list[0])
    for x in p2_list[1:3]:
        param_list.append(int(x))
    e_vals = list(p2_list[3].split(','))
    param_list.append(e_vals)
    for x in p2_list[4:]:
        param_list.append(os.path.join(p2_list[0],x)) 
    return param_list


def txtm_tarlen(loci_list_path, org_list_path, geneseq_path, out_path):
    """Generates a txt file that can be appended to a HypPiper output. 
    Prints the lengths of the genes for each txtm. """
    loci_list = filelines_to_list(loci_list_path)
    org_list = filelines_to_list(org_list_path)
    lens_dict = {}
    for org in org_list:
        org_lens_list = []
        for loci in loci_list:
            loci_rec_dict = SeqIO.to_dict(SeqIO.parse(
                f"{geneseq_path}/{loci}.fa", "fasta"))
            try:
                org_lens_list.append(len(loci_rec_dict[f"{org}-{loci}"].seq))
            except:
                org_lens_list.append('0')
            lens_dict[org] = org_lens_list
    with open(out_path, 'w') as out_handle:
        for item in lens_dict:
            line_string = ('\t').join(str(x) for x in lens_dict[item])
            out_handle.write(f"{item}\t{line_string}\n")
    return


def main():
    # Run the pipeline
    # args = sys.argv[1:]
    args = ["--param", "/Users/josec/Desktop/git_repos/TxtmFishing/Feb5-param.txt"]
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
    # for x in param_list:
    #     print(x)
    txtm_fishing_pipe(param_list)
 


if __name__ == "__main__":
    main()
