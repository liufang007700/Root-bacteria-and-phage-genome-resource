## CRBD BGC analysis

### 1. calculate contig length for each genome and subset BGC to contigs >=5Kb

#### a. Calculate contig length

```
seqkit fx2tab -l {1}.fna "|" cut -f 1,4 ">" summarise_genome_contig_length/{1}_length.txt \
 ::: $(cut -f 1 all_14242_genomes_ID.txt)
for i in  $(cut -f 1 all_14242_genomes_ID.txt)
do
sed -i "s/^/${i}\t/;s/ //g" summarise_genome_contig_length/${i}_length.txt
done
cat summarise_genome_contig_length/* > all_14242_genomes_contig_length.txt 
```

#### b. Summarize BGC from AntiSmash output

```

```

#### c. Edit BGC contigID with that of genome contigID
```
Please refer to the script file : 
```

#### c. Subset BGC to contigs >=5kb
```
cat <(awk '$4*1>=5000' all_14242_genomes_contig_length_up_partial_match_subset_up.txt| cut -f 1) <(awk '$3*1>=5000' all_14242_genomes_contig_length_up_exact_word_match_subset.txt | cut -f 2) > Iso_MAG_Pub_5kbplus_contigID_list
grep -w -F -f Iso_MAG_Pub_5kbplus_contigID_list Iso_MAG_Pub_BGC_cluster_info_updated.txt  > Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt
```

#### d. Subset BGC to CRBD NRgenome subset

```
grep -w -F -f CRBD_NRgenomeID_list  Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt > CRBD_NRgenome_5Kplus_BGC_cluster_info.txt
wc -l CRBD_NRgenomeID_list # 7531
cut -f 1 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt| sort | uniq | wc -l # 7278 has BGC within contig length >=5K
##### after double check, some of the genomes has BGC but their corresponding contig length is smaller than 5K
cut -f 2 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | sort | uniq | wc -l #55285 unique BGC together
cut -f 1,5 /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_bigscape/CRBD_NRgenome_5Kplus_bigscape_out/network_files/2023-09-06_23-15-07_auto_NRgenome_BGC_5kbplus/Network_Annotations_Full.tsv | grep -v '^BGC' > CRBD_NRgenome_5Kplus_BGC_bigscapeClass.txt # Pub_GCA_004305175.1_ASM430517v1_genomic  SIML01000028.1.region001.gbk    SIML01000028.1RhizobiumleguminosarumstrainSM86frag_SM86_22, 0   5717    RiPP    proteusin this gbk is included within the bigscape gbk input folder. however, within the network file this gbk is lost. So, I will manually add this line to the BGC_bigscapeClass file
cat CRBD_NRgenome_5Kplus_BGC_bigscapeClass.txt <(cat <(echo "SIML01000028.1.region001" "RiPP" | sed 's, ,\t,g')) > CRBD_NRgenome_5Kplus_BGC_bigscapeClass_up.txt
cut -f 2 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | sort | uniq | wc -l # 55285
wc -l CRBD_NRgenome_5Kplus_BGC_bigscapeClass_up.txt # 55285
cut -f 1-5 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | sort | uniq > CRBD_NRgenome_5Kplus_BGC_cluster_info_up.txt
## -------- now update the CRBD NRgenome BigscapeClass --------
### NOTE!!!! please be noted that within this CRBD_NRgenome_5Kplus_BGC_cluster_info_up.txt, there are three extra genomes and their BGCs, "Iso_Wt_W139-2"  "Iso_Wt_WTC65-1" "Iso_Wt_WTC65-2" 
```

#### e. calculate novel BGC versus bigScape

```
## prepare files

mkdir BGC_bigscape
cd BGC_bigscape
ln -s ../BGC_NRgenome/CRBD_NRgenome_5Kplus_BGC_cluster_info.txt ./                                                  
cut -f 1,2 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | sort | uniq > CRBD_NRgenomeID_gbk_map.txt
grep -w -F -f <(cut -f 1 6109_CRBD_HQ_NRgenome_meta.txt | grep -v 'GenomeID')  CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | cut -f 1,2 | grep -v 'Iso_Wt_W139-2' | grep -v 'Iso_Wt_WTC65-1' | grep -v 'Iso_Wt_WTC65-2' | sed 's,Iso_Wt_WTC89_2,Iso_Wt_WTC89-2,g'  >  CRBD_HQ_NRgenomeID_gbk_map.txt
grep -w -F -f  <(cut -f 1 CRBD_HQ_NRgenomeID_gbk_map.txt | sort | uniq)  <(cut -f 1,19-25 6109_CRBD_HQ_NRgenome_meta.txt | sed 's,d__,,g'| sed 's,p__,,g' |sed 's,c__,,g' | sed 's,o__,,g' | sed 's,f__,,g' | sed 's,g__,,g' | sed 's,s__,,g') > CRBD_HQ_NRgenome_5Kplus_BGC_harboring_taxonomy.txt

cd NRgenome_5kbplus_gbk
IFS=$'\n'
for file in $(cat ../CRBD_NRgenomeID_gbk_map.txt)
do
    echo $file
    genomeID=$(echo $file |cut -f 1)
    gbk=$(echo $file | cut -f 2)
    cp /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/"$genomeID"/"$gbk" ./
done

ls > temp
grep -v -w -F -f NRgenome_5kbplus_gbk/temp  <(cut -f 2 CRBD_NRgenomeID_gbk_map.txt)

## run bigScape

python ~/software/BIG_SCAPE/BiG-SCAPE-1.1.5/bigscape.py -l NRgenome_BGC_5kbplus -i /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_bigscape/NRgenome_5kbplus_gbk -o CRBD_NRgenome_5Kplus_bigscape_out_v3 --pfam_dir /mnt/m1/liufang/software/BIG_SCAPE/BiG-SCAPE-1.1.5/ -c 96 --include_gbk_str region --include_singletons --cutoffs 1.0 --clans-off --hybrids-off --mode auto --mibig --verbose >> run_v3.log 2>&1

python ~/software/BIG_SCAPE/BiG-SCAPE-1.1.5/bigscape.py -l HQ_NRgenome_BGC_5kbplus -i /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_bigscape/HQ_NRgenome_5kbplus_gbk -o CRBD_HQ_NRgenome_5Kplus_bigscape_out_no_MIBIG --pfam_dir /mnt/m1/liufang/software/BIG_SCAPE/BiG-SCAPE-1.1.5/ -c 96 --include_gbk_str region --include_singletons --cutoffs 0.65 --clans-off --hybrids-off --mode auto --verbose >> run_v4.log 2>&1
```
#### f. BGC BIGSCLICE

```
## Download BIG-FAM database

wget https://www.bioinformatics.nl/~kauts001/ltr/bigslice/paper_data/

## prepare the full input data for CRBD_HQ_NRgenome

cd /mnt/m4/liufang/db/BIG_SLICE
mkdir -p full_input_data/CRBD_HQ_NRgenome/
cd full_input_data/CRBD_HQ_NRgenome/
ln -s /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_bigscape/CRBD_HQ_NRgenomeID_gbk_map.txt ./
IFS=$'\n'

for file in $(cat CRBD_HQ_NRgenomeID_gbk_map.txt| cut -f 1| sort | uniq); do echo $file; mkdir -p full_input_data/CRBD_HQ_NRgenome/$file; done
for file in $(cat CRBD_HQ_NRgenomeID_gbk_map.txt); do genomeID=$(echo $file | cut -f 1); gbk=$(echo $file | cut -f 2); cp /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/"$genomeID"/"$gbk" full_input_data/CRBD_HQ_NRgenome/"$genomeID"; done

cp /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_14_length_127376_cov_66.943613.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/
cp /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_1_length_401843_cov_69.689376.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/
cp /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_9_length_176810_cov_67.604567.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/

bigslice -i /mnt/m4/liufang/db/BIG_SLICE/full_input_data /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output --threshold 900  -t 64 --program_db_folder /mnt/m1/liufang/anaconda3/envs/bigslice/bin/bigslice-models # the output are store in the ``CRBD_HQ_NRgenome_bigslice_output`` file
```

#### g. calculate the cos distance

```

```




