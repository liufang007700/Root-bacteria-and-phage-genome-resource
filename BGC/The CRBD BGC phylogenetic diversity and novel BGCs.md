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
#### b. Subset BGC to contigs >=5kb

```
Iso_MAG_Pub_BGC_cluster_info_updated.txt
awk '$3*1>=5000' all_14242_genomes_contig_length.txt > Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt
wc -l Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt  # 118481 BGCs
cut -f 1 Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt  | sort | uniq | wc -l # 13392 genomes

```

#### calculate novel BGC versus bigScape

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

```



