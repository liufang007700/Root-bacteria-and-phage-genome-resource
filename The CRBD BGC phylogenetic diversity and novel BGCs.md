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
#### b. Subset BGC to contigs >=5kb

```
awk '$3*1>=5000' all_14242_genomes_contig_length.txt > Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt
wc -l Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt  # 118481 BGCs
cut -f 1 Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt  | sort | uniq | wc -l # 13392 genomes
```

#### calculate novel BGC versus bigScape

```

```

