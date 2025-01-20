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
cd CRBD_HQ_NRgenome_bigslice_output
mkdir calculate_cosine_dist
### please be noted that BigFam BGC database include MiBIG2.0 (1910 BGC)

python convert_bigslice_data_2_table.py --input result/data.db -e False -c True
python extract_feature_matrix.py  ./ CRBD_HQ_NRgenome_BGC_feature_matrix
sed -i  -E 's,^\t,bgc_id\t,g' CRBD_HQ_NRgenome_BGC_feature_matrix
mv CRBD_HQ_NRgenome_BGC_feature_matrix  CRBD_plus_BiGfam_BGC_feature_matrix.txt

## Extract CRBD BGC feature matrix
csvtk filter -f 'dataset_id=11' /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > CRBD_bgc_id_list
#csvtk -t grep -P CRBD_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt  > sub_CRBD_BGC_feature_matrix.txt
csvtk -t grep -P CRBD_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt -o sub_CRBD_BGC_feature_matrix.txt -T

## Extract MiBIG BGC feature matrix

csvtk filter -f 'dataset_id=1' /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > MiBIG_bgc_id_list
csvtk -t grep -P MiBIG_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt -o sub_MiBIG_BGC_feature_matrix.txt -T

## Extract Other BGC feature (except CRBD) matrix

csvtk filter -f 'dataset_id!=11' /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > Other_bgc_id_list
csvtk -t grep -P Other_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt -o sub_Other_BGC_feature_matrix.txt -T

## calculate GCF membership

python convert_bigslice_data_2_table.py --input result/data.db -e False -c True
python extract_feature_matrix.py  ./ CRBD_HQ_NRgenome_BGC_feature_matrix
python clustering_bgc_from_bigslice_pool3.py sub_CRBD_BGC_feature_matrix.txt --prefix CRBD_cos_dist_based_clustering --maxBGCid 0
cat <(cut -f 1 sub_MiBIG_BGC_feature_matrix.txt | paste -s -d '\t')   <(paste -d '\t' <(cut -f 1 sub_CRBD_BGC_feature_matrix.txt) <(cut -f 2- CRBD_vs_MiBIG_cosine_distances_matrix.txt) | tail -n +2)  > CRBD_vs_MiBIG_cosine_distances_matrix_up.txt

python clustering_bgc_from_bigslice_pool3.py --input sub_CRBD_BGC_feature_matrix.txt --prefix CRB
D_self --maxBGCid 0
cut -f 1,2 CRBD_selfclustering.tsv | sed -E 's,\t,\tGCF_,g' | sed 's,GCF_cluster_0.2,GCF_id,g' > CRBD_HQ_NRgenome_BGC_GCF_membership.txt
```

#### h. Calculate the minimum cos distance between CRBD HQ NRgenome versus BigFam

```
cd /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files
CRBD_vs_BigFam_split_*.txt > /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/CRBD_vs_MiBiG_cos_dist/file_list
for file in $(cat file_list)
do
	echo $file
	ln -s /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/"$file" ./
done
parallel -j 100 /usr/bin/Rscript Summary_min_cos_dist.R {1}_cosine_distances_matrix.txt ./min_cos_dist/{1} ::: $(sed 's,_cosine_distances_matrix.txt,,g' file_list)
grep -v -w -F -f file_list  <(ls /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/CRBD_vs_BigFam_split_*.txt | sed 's,/mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/,,g') > left_file_list
for file in $(cat left_file_list)
do
    echo $file
    ln -s /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/"$file" ./
done
parallel -j 100 /usr/bin/Rscript Summary_min_cos_dist.R {1}_cosine_distances_matrix.txt ./min_cos_dist/{1} ::: $(sed 's,_cosine_distances_matrix.txt,,g' left_file_list)

ulimit -n 4096 # other wide paste command will report error like "too many files opened"
paste -d "\t" min_cos_dist/* > combine_split_min_cos_dist.txt
paste -d "\t" <(cut -f 1 combine_split_min_cos_dist.txt) <(cat  <(ls min_cos_dist/* | sed 's,min_cos_dist/CRBD_vs_BigFam_,,g' | sed 's,_min_cos_dist.txt,,g' | paste -d '\t' -s)  <(cut -f $(seq 2 2 $(head -2 combine_split_min_cos_dist.txt | tail -1 | awk '{print NF}') | paste -s -d ',') combine_split_min_cos_dist.txt | tail -n +2 | sed 's, ,\t,g')) | sed 's,Min_cosdist,bgc_id,g' > CRBD_vs_BigFam_min_cos_dist_combine.txt
mkdir split_cosine_distances_matrix
mv CRBD_vs_BigFam_split_*_cosine_distances_matrix.txt split_cosine_distances_matrix

cut -d ',' -f 2,9 bgc.csv > BGC_name_vs_bgc_id_map.txt

## Summary_min_cos_dist.R

#!/bin/Rscript

library(stringr)
library(reshape2)
library(dplyr)

#Reading parameters
infile <- commandArgs(trailingOnly=T)[1];
outfile_prefix <- commandArgs(trailingOnly=T)[2];
dist_df<-read.table(infile,sep="\t",header=TRUE,row.names=1,quote="")
dist_df_up<-dist_df[apply(dist_df,1,sum)>0,apply(dist_df,2,sum)>0]
min_per_row_up<-data.frame(Min_cosdist=apply(dist_df_up,1,function(x) min(x)))
write.table(min_per_row_up,file=paste(outfile_prefix,"_min_cos_dist.txt",sep=""),sep="\t",quote=FALSE)
```

#### i. Generate the novel GCF composition and diversity

```
#!/bin/Rscript

#### after individual CRBD vs BigFam split min cos dist combined into one file, now calculate the minimum distance bewteen CRBD vs BigFam

library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)


taxa_color<-read.csv("GTDB_phylum_color_Ray_define_top12_Sep_05_2023.csv",header = FALSE,col.names = c("PhyClass","Color"))
rownames(taxa_color)<-taxa_color$PhyClass
PhyClass_order<-c("Alphaproteobacteria","Gammaproteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Myxococcota","Patescibacteria","Spirochaetota","Fibrobacterota","Chloroflexota","Verrucomicrobiota","Acidobacteriota","Others")
BGC_color<-c("#e59572","#65a9b8","#d9a7a7","#bf7e7e","#a86060","#d4cb7b","#d1a7db","#5fa171","grey")
BGC_class<-c("NRPS","Others","PKS-NRP_Hybrids","PKSI","PKSother","RiPPs","Saccharides","Terpene","Multiple_BigScapeClasses")
BGC_class_annot<-data.frame(Color=BGC_color,Class=BGC_class)


## calculate the min cos distance between CRBD BGC against all MiBIG BGC

dist_df<-read.table("CRBD_vs_BigFam_min_cos_dist_combine.txt",sep="\t",header=TRUE,row.names=1,quote="")
dim(dist_df) # 48635  1226
min_per_row_up<-data.frame(Min_cosdist=apply(dist_df,1,function(x) min(x)))

## now read into the CRBD GCF member dataframe
GCF_mem<-read.table("CRBD_HQ_NRgenome_BGC_GCF_membership.txt",header=TRUE,sep="\t")
CRBD_GCF_mem<-GCF_mem%>%filter(bgc_id%in%row.names(min_per_row_up)) # 48635
CRBD_GCF_mem$bgc_id<-as.character(CRBD_GCF_mem$bgc_id)
CRBD_GCF_mem$gcf_id<-as.character(CRBD_GCF_mem$gcf_id)
CRBD_vs_BigFam_min_cosdist_GCF_mem<-inner_join(data.frame(bgc_id=row.names(min_per_row_up),min_per_row_up),CRBD_GCF_mem%>%select(c(bgc_id,gcf_id)))
length(unique(CRBD_vs_BigFam_min_cosdist_GCF_mem$gcf_id)) # 12865 GCF
CRBD_vs_BigFam_min_cosdist_GCF_mean<-data.frame(CRBD_vs_BigFam_min_cosdist_GCF_mem%>%group_by(gcf_id)%>%mutate(min_cosdist_GCF_mean=mean(Min_cosdist))%>%select(gcf_id,min_cosdist_GCF_mean)%>%unique())
summary(CRBD_vs_BigFam_min_cosdist_GCF_mean$min_cosdist_GCF_mean) # Min 1st Median Mean 3rd Max 0.01746 0.30399 0.45248 0.45536 0.60384 0.86563 

## GCF BigScapeClass composition --- cutoff scheme 1

range_cutoffs <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)
CRBD_vs_BigFam_min_cosdist_GCF_mean$range_category<-cut(CRBD_vs_BigFam_min_cosdist_GCF_mean$min_cosdist_GCF_mean,breaks = range_cutoffs, include.lowest = TRUE)
write.table(CRBD_vs_BigFam_min_cosdist_GCF_mean,"R_write_CRBD_vs_BigFam_mean_cos_dist_perGCF.txt",sep="\t",row.names=FALSE,quote=FALSE) 
## prepare gcf_id and BGC bigscape class mapping file

bgc_id_vs_name<-read.csv("bgc.csv",header=TRUE)
BigSlice_bgc_id_vs_name<-bgc_id_vs_name%>%select(c(id,orig_filename))%>%dplyr::rename(bgc_id=id)
BigSlice_bgc_id_vs_name$bgc_id<-factor(BigSlice_bgc_id_vs_name$bgc_id)
BigSlice_bgc_id_vs_name_add_gcf_id<-inner_join(BigSlice_bgc_id_vs_name,CRBD_GCF_mem%>%select(gcf_id,bgc_id))
dim(BigSlice_bgc_id_vs_name_add_gcf_id) # 48635 since 9 of the BGC does not has cos distance data
BigSlice_bgc_id_vs_name_add_gcf_id$BGC<-gsub(".gbk","",BigSlice_bgc_id_vs_name_add_gcf_id$orig_filename)
CRBD_BGC_bigscape_class<-read.table("CRBD_HQ_NRgenome_BGC_bigscape_class.txt",sep="\t",header=TRUE)
BigSlice_BGC_GCF_bigscape_class<-inner_join(BigSlice_bgc_id_vs_name_add_gcf_id,CRBD_BGC_bigscape_class)
BigSlice_BGC_GCF_class_summary<-data.frame(BigSlice_BGC_GCF_bigscape_class%>%group_by(gcf_id,BiG.SCAPE.class)%>%mutate(BGC_SUM_perGCF_perBigScapeClass=n())%>%select(c(gcf_id,BiG.SCAPE.class,BGC_SUM_perGCF_perBigScapeClass))%>%unique())# It turned out that some of the GCF belongs to different BIG_scape classes, to solve this problem, the BigClass of the GCF would be assigned by the bigscape class habors the maximum number of BGCs
BigSlice_BGC_GCF_BigScapeClass<-BigSlice_BGC_GCF_class_summary%>%group_by(gcf_id)%>%mutate(GCF_max_BigScapeClass=max(BGC_SUM_perGCF_perBigScapeClass))%>%filter(GCF_max_BigScapeClass==BGC_SUM_perGCF_perBigScapeClass) # 6324 rows, while unique GCF are 6276, for this reason, I decided to assign multiple_classes to GCF assigned to several different bigscape classes
BigSlice_BGC_GCF_class_summary_up<-BigSlice_BGC_GCF_class_summary%>%group_by(gcf_id)%>%mutate(BigScapeClass_SUM_per_GCF=n())%>%mutate(BigScapeClass_up=if_else(BigScapeClass_SUM_per_GCF>1,"Multiple_BigScapeClasses",BiG.SCAPE.class))%>%select(c(gcf_id,BigScapeClass_up))%>%unique() # 12865 GCFs

## now the BigScape Class for each GCF is ready too

CRBD_vs_BigFam_min_cosdist_GCF_mean_BigScapeClass<-inner_join(CRBD_vs_BigFam_min_cosdist_GCF_mean,BigSlice_BGC_GCF_class_summary_up)
CRBD_vs_BigFam_min_cosdist_GCF_mean_BigScapeClass$Value=1
CRBD_vs_BigFam_cosdist_distribution_BigScapeClass<-data.frame(CRBD_vs_BigFam_min_cosdist_GCF_mean_BigScapeClass%>%group_by(range_category,BigScapeClass_up)%>%mutate(Sum_perRange=sum(Value))%>%select(c(range_category,BigScapeClass_up,Sum_perRange))%>%unique())
CRBD_vs_BigFam_cosdist_distribution_BigScapeClass$BigScapeClass_up<-factor(CRBD_vs_BigFam_cosdist_distribution_BigScapeClass$BigScapeClass_up,levels=BGC_class) 

pdf("BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_cosdist_distribution_BigScapeClass,aes(x=range_category,y=Sum_perRange,fill=BigScapeClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=BGC_color)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()


pdf("Novel_GCF_only_BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_cosdist_distribution_BigScapeClass%>%filter(!range_category%in%c("[0,0.05]","(0.05,0.1]","(0.1,0.15]","(0.15,0.2]")),aes(x=range_category,y=Sum_perRange,fill=BigScapeClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=BGC_color)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()

pdf("Except_1range_BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_cosdist_distribution_BigScapeClass%>%filter(!range_category%in%c("[0,0.05]")),aes(x=range_category,y=Sum_perRange,fill=BigScapeClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=BGC_color)+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()


############ GCF range  scheme 1 ################

## GCF taxa summary

BGC_taxa<-read.table("R_write_CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa.txt",sep="\t",header=TRUE,quote="")
dim(BigSlice_bgc_id_vs_name_add_gcf_id)
BigSlice_bgc_id_vs_name_add_gcf_id[1:5,]
BigSlice_bgc_id_vs_name_add_gcf_id_taxa<-inner_join(BigSlice_bgc_id_vs_name_add_gcf_id%>%select(!BGC)%>%dplyr::rename(BGC=orig_filename),BGC_taxa%>%select(c(BGC,PhyClass_collapse)))
GCF_PhyClass_SUM<-data.frame(BigSlice_bgc_id_vs_name_add_gcf_id_taxa%>%group_by(gcf_id,PhyClass_collapse)%>%mutate(GCF_PhyClass_SUM=n())%>%select(c(gcf_id,PhyClass_collapse,GCF_PhyClass_SUM))%>%unique())
GCF_PhyClass<-GCF_PhyClass_SUM%>%group_by(gcf_id)%>%mutate(GCF_PhyClass_num=n())%>%select(!(GCF_PhyClass_SUM))%>%unique()%>%mutate(PhyClass_up=if_else(GCF_PhyClass_num>1,"Multiple_PhyClass",PhyClass_collapse))%>%select(!PhyClass_collapse)%>%unique() # 6276 GCFs

## combine GCF taxa with mean cos dist per GCF
range_cutoffs <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)
CRBD_vs_BigFam_min_cosdist_GCF_mean$range_category<-cut(CRBD_vs_BigFam_min_cosdist_GCF_mean$min_cosdist_GCF_mean,breaks = range_cutoffs, include.lowest = TRUE)
CRBD_vs_BigFam_min_cosdist_GCF_mean_add_taxa<-inner_join(CRBD_vs_BigFam_min_cosdist_GCF_mean,GCF_PhyClass)
CRBD_vs_BigFam_GCF_taxa_summary<-data.frame(CRBD_vs_BigFam_min_cosdist_GCF_mean_add_taxa%>%group_by(range_category,PhyClass_up)%>%mutate(PhyClassSum_perRange=n())%>%select(range_category,PhyClassSum_perRange,PhyClass_up)%>%unique())
CRBD_vs_BigFam_GCF_taxa_summary$PhyClass_up<-factor(CRBD_vs_BigFam_GCF_taxa_summary$PhyClass_up,levels=c(PhyClass_order[-7],"Multiple_PhyClass"))

pdf("BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges_taxonomy.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_GCF_taxa_summary,aes(x=range_category,y=PhyClassSum_perRange,fill=PhyClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=c(taxa_color[PhyClass_order[-7],]$Color,"#5e5e5e"))+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()

pdf("Novel_GCF_only_BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges_taxonomy.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_GCF_taxa_summary%>%filter(!range_category%in%c("[0,0.05]","(0.05,0.1]","(0.1,0.15]","(0.15,0.2]")),aes(x=range_category,y=PhyClassSum_perRange,fill=PhyClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=c(taxa_color[PhyClass_order[-7],]$Color,"#5e5e5e"))+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()

pdf("Except_1range_BigSlice_GCF_CRBD_HQ_NRgenome_GCF_vs_BigFam_cosine_distance_along_ranges_taxonomy.pdf",width=8,height=4)
ggplot(CRBD_vs_BigFam_GCF_taxa_summary%>%filter(!range_category%in%c("[0,0.05]")),aes(x=range_category,y=PhyClassSum_perRange,fill=PhyClass_up))+geom_bar(stat = "identity",size=0.5,color="#f0f2f2",width = 0.8)+scale_fill_manual(values=c(taxa_color[PhyClass_order[-7],]$Color,"#5e5e5e"))+theme_bw()+theme(axis.title.x = element_blank(),axis.text.x=element_text(angle=90),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+labs(y="GCF number",x="Distance cutoff")
dev.off()
```






