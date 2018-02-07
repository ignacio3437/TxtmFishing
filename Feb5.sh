########Pulling out all of the genes that were targeted from transcriptome data



#Set proj_home
export proj_home="/Users/josec/Desktop/Crinoid_capture/Feb5_hybTxCrinoid"

# python3 /Users/josec/Desktop/git_repos/TxtmFishing/TxtmFishing.py --param Feb5_Strict_tarlens.txt


# cd $proj_home/Exonerate_StrictFeb5_clean
# for x in *.fa;do cat $x >> StrictFeb5_targets.fa;done
# mv StrictFeb5_targets.fa ..


#Run HybPiper
cd $proj_home
mkdir HybPiper_Strict
cp Cri_names.txt HybPiper_Strict
cp StrictFeb5_targets.fa HybPiper_Strict
cd $proj_home/HybPiper_Strict
while read name
do /Users/josec/Desktop/Gitclones/HybPiper/reads_first.py -r ../RawReads/$name*.fastq -b StrictFeb5_targets.fa --prefix $name --bwa
done < Cri_names.txt

#Get the lengths of contigs recovered.
python /Users/josec/Desktop/Gitclones/HybPiper/get_seq_lengths.py StrictFeb5_targets.fa Cri_names.txt dna > $proj_home/StrictFeb5_seq_lengths.txt
#Generate Summary Data
python /Users/josec/Desktop/Gitclones/HybPiper/hybpiper_stats.py $proj_home/StrictFeb5_seq_lengths.txt Cri_names.txt > $proj_home/StrictFeb5_stats.txt
#Pull DNA seqs
mkdir $proj_home/StrictFeb5_dna_seqs
cd $proj_home/StrictFeb5_dna_seqs
python /Users/josec/Desktop/Gitclones/HybPiper/retrieve_sequences.py $proj_home/StrictFeb5_targets.fa $proj_home/HybPiper_Strict dna

#Get seq lengths of txtm terminals
cd /Users/josec/Desktop/git_repos/TxtmFishing
mkdir $proj_home/StrictFeb5_bygene
python3 -c "import TxtmFishing; TxtmFishing.dir_CATer('$proj_home/Exonerate_clean','$proj_home/StrictFeb5_dna_seqs','$proj_home/StrictFeb5_bygene')"
python3 -c "import TxtmFishing; TxtmFishing.txtm_tarlen('$proj_home/Loc_list2.txt','$proj_home/TranscriptomeList.txt','$proj_home/Exonerate_StrictFeb5_clean','$proj_home/StrictFeb5_Txm_tarlens.txt')"

#append this to end of StrictFeb5_seq_lengths.txt file 




################################
# cd $proj_home/Exonerate_LooseFeb5_clean
# for x in *.fa;do cat $x >> LooseFeb5_targets.fa;done
# mv LooseFeb5_targets.fa ..

#Run HybPiper
cd $proj_home
mkdir HybPiper_Loose
cp Cri_names.txt HybPiper_Loose
cp LooseFeb5_targets.fa HybPiper_Loose
cd $proj_home/HybPiper_Loose
while read name
do /Users/josec/Desktop/Gitclones/HybPiper/reads_first.py -r ../RawReads/$name*.fastq -b LooseFeb5_targets.fa --prefix $name --bwa
done < Cri_names.txt

#Get the lengths of contigs recovered.
python /Users/josec/Desktop/Gitclones/HybPiper/get_seq_lengths.py LooseFeb5_targets.fa Cri_names.txt dna > $proj_home/LooseFeb5_seq_lengths.txt
#Generate Summary Data
python /Users/josec/Desktop/Gitclones/HybPiper/hybpiper_stats.py $proj_home/LooseFeb5_seq_lengths.txt Cri_names.txt > $proj_home/LooseFeb5_stats.txt
#Pull DNA seqs
mkdir $proj_home/LooseFeb5_dna_seqs
cd $proj_home/LooseFeb5_dna_seqs
python /Users/josec/Desktop/Gitclones/HybPiper/retrieve_sequences.py $proj_home/LooseFeb5_targets.fa $proj_home/HybPiper_Loose dna

#Get seq lengths of txtm terminals
cd /Users/josec/Desktop/git_repos/TxtmFishing
mkdir $proj_home/LooseFeb5_bygene
python3 -c "import TxtmFishing; TxtmFishing.dir_CATer('$proj_home/Exonerate_clean','$proj_home/LooseFeb5_dna_seqs','$proj_home/LooseFeb5_bygene')"
python3 -c "import TxtmFishing; TxtmFishing.txtm_tarlen('$proj_home/Loc_list2.txt','$proj_home/TranscriptomeList.txt','$proj_home/Exonerate_LooseFeb5_clean','$proj_home/LooseFeb5_Txm_tarlens.txt')"

#append this to end of LooseFeb5_seq_lengths.txt file 


################################



# cd $proj_home/HybTxt_seqs/aln
# # for fna in *.FNA;do iqtree -s $fna -nt AUTO -m TEST;done
# parallel --jobs 4 iqtree -s {} -nt 2 -m TEST ::: *.FNA
# cat *.treefile > /Users/josec/Desktop/Crinoid_capture/Feb_HybTxCrinoid/HybTxt_Astral/Feb_genetrees.tre




# #ASTRAL
# cd $proj_home/Hybtxt_Astral
# java -jar /Users/josec/Desktop/Gitclones/ASTRAL/astral.5.5.6.jar -i genetrees.tre -o Feb2_astral 2>Feb2_astral.log





