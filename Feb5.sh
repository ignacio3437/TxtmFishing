########Pulling out all of the genes that were targeted from transcriptome data



#Set proj_home
export proj_home="/Users/josec/Desktop/Crinoid_capture/Feb5_hybTxCrinoid"

python3 /Users/josec/Desktop/git_repos/TxtmFishing/TxtmFishing.py --param Feb5-param.txt


cd $proj_home/Txtms_TrimtoTargets_alignments
for x in *.fa;do cat $x >> Feb2_Trim_targets.fa;done
mv Feb2_Trim_Targets.fa ..


#Run HybPiper
cd $proj_home
mkdir HybPiper
cp Strict_cri_names.txt HybPiper
cp Feb2_Trim_targets.fa HybPiper
cd $proj_home/HybPiper
while read name
do /Users/josec/Desktop/Gitclones/HybPiper/reads_first.py -r ../RawReads/$name*.fastq -b Feb2_Trim_targets.fa --prefix $name --bwa
done < Strict_cri_names.txt
# done < Strict_cri_names2.txt

#Get the lengths of contigs recovered.
python /Users/josec/Desktop/Gitclones/HybPiper/get_seq_lengths.py Feb2_Trim_targets.fa Strict_cri_names.txt dna > $proj_home/test_seq_lengths.txt
#Generate Summary Data
python /Users/josec/Desktop/Gitclones/HybPiper/hybpiper_stats.py $proj_home/test_seq_lengths.txt Strict_cri_names.txt > $proj_home/test_stats.txt

mkdir $proj_home/dna_seqs
cd $proj_home/dna_seqs
python /Users/josec/Desktop/Gitclones/HybPiper/retrieve_sequences.py $proj_home/Feb2_Trim_targets.fa $proj_home/HybPiper dna


cd $proj_home/bin
mkdir $proj_home/Feb2_bygene
python3 -c "import TxtmFishing; TxtmFishing.dir_CATer('$proj_home/Exonerate_clean','$proj_home/dna_seqs','$proj_home/Feb2_bygene')"


python3 -c "import TxtmFishing; TxtmFishing.txtm_tarlen('$proj_home/Loc_list2.txt','$proj_home/TranscriptomeList.txt','$proj_home/Exonerate_clean','$proj_home/Txm_tarlens.txt')"
# python3 -c "import TxtmFishing; TxtmFishing.txtm_tarlen('$proj_home/Trim_Loc_List.txt','$proj_home/TranscriptomeList.txt','$proj_home/Exonerate_clean','$proj_home/Txm_tarlens.txt')"

#append this to end of test_seq_lengths.txt file 








################################


cd $proj_home/bin
for fasta_file in $proj_home/ExonerateOut/*.fa; do python3 -c "import Feb2; Feb2.fasta_deduper('$fasta_file')";done

$proj_dir/Update2/Gene_aln


#align with mafft
cd $proj_home/HybTxt_seqs/
mkdir aln
for fna in *.FNA; do mafft --thread 8 --auto $fna > $proj_home/HybTxt_seqs/aln/$fna;done

#Trim alingments
mkdir aln_trim
cd aln
for fna in *.FNA; do trimal -strict -in $fna -out ../aln_trim/$fna;done


mkdir $proj_home/HybTxt_concat
./catfasta2phyml.pl -c -f $proj_home/Update2/Gene_aln/*.FNA > $proj_home/HybTxt_concat/Feb2_cat.fasta
trimal -strict -in Feb2_cat.fasta -out TRIM_Feb2_cat.fasta



cd $proj_home/HybTxt_seqs/aln
# for fna in *.FNA;do iqtree -s $fna -nt AUTO -m TEST;done
parallel --jobs 4 iqtree -s {} -nt 2 -m TEST ::: *.FNA
cat *.treefile > /Users/josec/Desktop/Crinoid_capture/Feb_HybTxCrinoid/HybTxt_Astral/Feb_genetrees.tre



#trim gene trees
cd $proj_home/HybTxt_seqs/aln
parallel --jobs 4 iqtree -s {}  -nt 2 -m TEST ::: *.FNA
cat *.treefile > /Users/josec/Desktop/Crinoid_capture/Feb_HybTxCrinoid/HybTxt_Astral/Feb_TRIM_genetrees.tre




#ASTRAL
cd $proj_home/Hybtxt_Astral
java -jar /Users/josec/Desktop/Gitclones/ASTRAL/astral.5.5.6.jar -i genetrees.tre -o Feb2_astral 2>Feb2_astral.log

echo 'done'







################################


mkdir $proj_home/dna_seqs
cd $proj_home/dna_seqs
python /Users/josec/Desktop/Gitclones/HybPiper/retrieve_sequences.py $proj_home/Feb2_targets.fa .. dna
mkdir aln
for fna in *.FNA; do mafft --thread 8 --auto $fna > $proj_home/dna_seqs/aln/$fna;done

python $proj_home/blastParser3.py
