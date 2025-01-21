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
$proj=/mnt/m2/dairui/project/binning
## combine BGC results
cd $proj/MAG_finalization/all
mkdir -p 14_BGC/sep_genome
cd $proj/MAG_finalization/all/14_BGC/sep_genome
ln -s $proj/isolate/08_BGCs/Arabidopsis/sep_genome/* ./
ln -s $proj/isolate/08_BGCs/Maize/sep_genome/* ./
ln -s $proj/isolate/08_BGCs/Medicago/sep_genome/* ./
ln -s $proj/isolate/08_BGCs/Rice/sep_genome/* ./
ln -s $proj/isolate/08_BGCs/Wheat/sep_genome/* ./
ln -s $proj/MAG_finalization/*/08_BGCs/sep_genome/* ./
ln -s $proj/MAG_finalization/Rc_HNmini/*/08_BGCs/sep_genome/* ./
ln -s $proj/published_genomes/Pub_3556_antismash/sep_genome/* ./
ln -s $proj/published_genomes/IMG_root_genomes_0330/08_antismash/sep_genome/* ./
cd $proj/MAG_finalization/all/14_BGC/
ls sep_genome > all_14242_genomes_ID

DIR=$proj/MAG_finalization/all/14_BGC/sep_genome
cd $DIR/..

for i in `cat $DIR/../all_14242_genomes_ID`; do
for j in `ls sep_genome/$i/*region*.gbk`;do
a=`ls $j | cut -f 3 -d '/'`
for n in `seq 1 $(grep -c "category" $j)`; do
echo -n -e "$i\t$a\t" >> cluster.list; \
grep -P 'DEFINITION|Orig\.' $j |sed 's/  *//g;s/DEFINITION//;s/\.$//;s/Orig.start:://;s/Orig.end:://' |tr '\n' '\t' >> cluster.list
grep "/" $j | grep -A 6 -m $n 'category'| tail -n 7| grep -E "category|product" | sed 's/                     \/category="//;s/                     \/product="//'|tr '\n' '\t'|cut -f 1,2 -d '"' | sed 's/"//'>> cluster.list
done
done
done &
sort cluster.list | uniq  > cluster_info.txt
```

#### c. Edit BGC contigID with that of genome contigID
```
Please refer to the script file : ``https://github.com/liufang007700/Root-bacteria-and-phage-genome-resource/blob/main/BGC/CODE/match_contigID_with_BGC_clusterID.sh``
```

#### c. Subset BGC to contigs >=5kb
```
cat <(awk '$4*1>=5000' all_14242_genomes_contig_length_up_partial_match_subset_up.txt| cut -f 1) <(awk '$3*1>=5000' all_14242_genomes_contig_length_up_exact_word_match_subset.txt | cut -f 2) > Iso_MAG_Pub_5kbplus_contigID_list
grep -w -F -f Iso_MAG_Pub_5kbplus_contigID_list Iso_MAG_Pub_BGC_cluster_info_updated.txt  > Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt
```

#### d. Subset BGC to CRBD NRgenome subset

```
wd=/mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA
grep -w -F -f CRBD_NRgenomeID_list  Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt > CRBD_NRgenome_5Kplus_BGC_cluster_info.txt
wc -l CRBD_NRgenomeID_list # 7531
cut -f 1 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt| sort | uniq | wc -l # 7278 has BGC within contig length >=5K
##### after double check, some of the genomes has BGC but their corresponding contig length is smaller than 5K
cut -f 2 CRBD_NRgenome_5Kplus_BGC_cluster_info.txt | sort | uniq | wc -l #55285 unique BGC together
cut -f 1,5 $wd/CRBD_BGC/BGC_bigscape/CRBD_NRgenome_5Kplus_bigscape_out/network_files/2023-09-06_23-15-07_auto_NRgenome_BGC_5kbplus/Network_Annotations_Full.tsv | grep -v '^BGC' > CRBD_NRgenome_5Kplus_BGC_bigscapeClass.txt # Pub_GCA_004305175.1_ASM430517v1_genomic  SIML01000028.1.region001.gbk    SIML01000028.1RhizobiumleguminosarumstrainSM86frag_SM86_22, 0   5717    RiPP    proteusin this gbk is included within the bigscape gbk input folder. however, within the network file this gbk is lost. So, I will manually add this line to the BGC_bigscapeClass file
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
    cp $proj/MAG_finalization/all/14_BGC/sep_genome/"$genomeID"/"$gbk" ./
done

ls > temp
grep -v -w -F -f NRgenome_5kbplus_gbk/temp  <(cut -f 2 CRBD_NRgenomeID_gbk_map.txt)

## run bigScape

python ~/software/BIG_SCAPE/BiG-SCAPE-1.1.5/bigscape.py -l NRgenome_BGC_5kbplus -i $wd/CRBD_BGC/BGC_bigscape/NRgenome_5kbplus_gbk -o CRBD_NRgenome_5Kplus_bigscape_out_v3 --pfam_dir /mnt/m1/liufang/software/BIG_SCAPE/BiG-SCAPE-1.1.5/ -c 96 --include_gbk_str region --include_singletons --cutoffs 1.0 --clans-off --hybrids-off --mode auto --mibig --verbose >> run_v3.log 2>&1

python ~/software/BIG_SCAPE/BiG-SCAPE-1.1.5/bigscape.py -l HQ_NRgenome_BGC_5kbplus -i $wd/CRBD_BGC/BGC_bigscape/HQ_NRgenome_5kbplus_gbk -o CRBD_HQ_NRgenome_5Kplus_bigscape_out_no_MIBIG --pfam_dir /mnt/m1/liufang/software/BIG_SCAPE/BiG-SCAPE-1.1.5/ -c 96 --include_gbk_str region --include_singletons --cutoffs 0.65 --clans-off --hybrids-off --mode auto --verbose >> run_v4.log 2>&1
```
#### f. BGC BIGSCLICE

```
## Download BIG-FAM database

wget https://www.bioinformatics.nl/~kauts001/ltr/bigslice/paper_data/

## prepare the full input data for CRBD_HQ_NRgenome

cd /mnt/m4/liufang/db/BIG_SLICE
mkdir -p full_input_data/CRBD_HQ_NRgenome/
cd full_input_data/CRBD_HQ_NRgenome/
ln -s $wd/CRBD_BGC/BGC_bigscape/CRBD_HQ_NRgenomeID_gbk_map.txt ./
IFS=$'\n'

for file in $(cat CRBD_HQ_NRgenomeID_gbk_map.txt| cut -f 1| sort | uniq); do echo $file; mkdir -p full_input_data/CRBD_HQ_NRgenome/$file; done
for file in $(cat CRBD_HQ_NRgenomeID_gbk_map.txt); do genomeID=$(echo $file | cut -f 1); gbk=$(echo $file | cut -f 2); cp $proj/MAG_finalization/all/14_BGC/sep_genome/"$genomeID"/"$gbk" full_input_data/CRBD_HQ_NRgenome/"$genomeID"; done

cp $proj/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_14_length_127376_cov_66.943613.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/
cp $proj/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_1_length_401843_cov_69.689376.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/
cp $proj/MAG_finalization/all/14_BGC/sep_genome/Iso_Wt_WTC89-2/Iso_Wt_WTC89-2___NODE_9_length_176810_cov_67.604567.region001.gbk full_input_data/CRBD_HQ_NRgenome/Iso_Wt_WTC89_2/

bigslice -i /mnt/m4/liufang/db/BIG_SLICE/full_input_data $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output --threshold 900  -t 64 --program_db_folder /mnt/m1/liufang/anaconda3/envs/bigslice/bin/bigslice-models # the output are store in the ``CRBD_HQ_NRgenome_bigslice_output`` file
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
csvtk filter -f 'dataset_id=11' $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > CRBD_bgc_id_list
#csvtk -t grep -P CRBD_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt  > sub_CRBD_BGC_feature_matrix.txt
csvtk -t grep -P CRBD_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt -o sub_CRBD_BGC_feature_matrix.txt -T

## Extract MiBIG BGC feature matrix

csvtk filter -f 'dataset_id=1' $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > MiBIG_bgc_id_list
csvtk -t grep -P MiBIG_bgc_id_list -f 1 CRBD_plus_BiGfam_BGC_feature_matrix.txt -o sub_MiBIG_BGC_feature_matrix.txt -T

## Extract Other BGC feature (except CRBD) matrix

csvtk filter -f 'dataset_id!=11' $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/result/Saved_Dataframes_data/bgc.csv | csvtk cut -f id | grep -v 'id' > Other_bgc_id_list
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
cd $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files
CRBD_vs_BigFam_split_*.txt > $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/CRBD_vs_MiBiG_cos_dist/file_list
for file in $(cat file_list)
do
	echo $file
	ln -s $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/"$file" ./
done
parallel -j 100 /usr/bin/Rscript Summary_min_cos_dist.R {1}_cosine_distances_matrix.txt ./min_cos_dist/{1} ::: $(sed 's,_cosine_distances_matrix.txt,,g' file_list)
grep -v -w -F -f file_list  <(ls $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/CRBD_vs_BigFam_split_*.txt | sed 's,$wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/,,g') > left_file_list
for file in $(cat left_file_list)
do
    echo $file
    ln -s $wd/CRBD_BGC/BGC_BiG_SliCE/CRBD_HQ_NRgenome_bigslice_output/calculate_cosine_dist/Other_BGC_split_1000_files/"$file" ./
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

## Rscript for main figure - Generate the novel GCF composition and diversity 

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
## Rscript for main figure  - BGC composition along the CRBD phylogenetic tree visualized at genus level

```
---
title: "Phylogenetic_distribution_of_BGC_cluster"
output: html_document
date: '2023-10-11'
editor_options: 
  chunk_output_type: console
---

# load envs and setup main_theme

library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(vegan)
library(phyloseq)
library(ggpubr)
library(pheatmap)
library(DESeq2)
library(directlabels)
library(tidyr)
library(extrafont)
library(circlize)
#font_import(pattern = "Arial") # this only needed once and no need to import every time we run loadfonts(device='pdf')
loadfonts(device = "pdf")
library(RColorBrewer)
library(tidyheatmap)
library(cowplot)
library(stringr)
library(patchwork)
main_theme<-theme(panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 15,face="bold"), axis.text.x = element_text(family="Arial",size=13,angle = 90,hjust = 0,vjust = 0),axis.text.y = element_text(family="Arial",size=13),axis.title.x = element_text(family="Arial",size=15,face="bold"),axis.title.y = element_text(family="Arial",size=15,face="bold"), legend.title = element_text(family="Arial",size=12,face="bold"),legend.text =element_text(family="Arial",size=12),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.box.background = element_rect(color="transparent", fill = "transparent"),legend.key = element_rect(fill = "transparent", colour = NA),strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),plot.margin=unit(c(1,1,1,1),"cm"))

taxa_color<-read.csv("/Users/fangliu/Documents/IGDB_Bai_lab/Database_and_color_pallete/Color_scheme/GTDB_phylum_color_Ray_define_top12_Sep_05_2023.csv",header = FALSE,col.names = c("PhyClass","Color"))
rownames(taxa_color)<-taxa_color$PhyClass

BGC_color<-c("#e59572","#65a9b8","#d9a7a7","#bf7e7e","#a86060","#d4cb7b","#d1a7db","#5fa171","grey")
BGC_class<-c("NRPS","Others","PKS-NRP_Hybrids","PKSI","PKSother","RiPPs","Saccharides","Terpene","Multiple_BigScapeClasses")
BGC_class_annot<-data.frame(Color=BGC_color,Class=BGC_class)

# Subset CRBD HQ NRgenome taxonomy

CRBD_genome_meta<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/9772_CRBD_genomes_metadata_0918.txt",sep="\t",header=TRUE,quote = "") #9772
summary(CRBD_genome_meta$GenomeID=="Pub_GCA_000468095.1_Pantoea_sp._AS-PWVM4_genomic") # TRUE=1
colnames(CRBD_genome_meta)<-gsub(pattern = "[.]",replacement = "_",colnames(CRBD_genome_meta))
colnames(CRBD_genome_meta)<-gsub(pattern = "__ncbi_","_ncbi",colnames(CRBD_genome_meta))
colnames(CRBD_genome_meta)<-gsub(pattern = "__gtdb_","_gtdb",colnames(CRBD_genome_meta))
CRBD_genome_meta_up<-inner_join(CRBD_genome_meta,CRBD_genome_meta%>%filter(Rep_NR_Genomes=="Yes")%>%select(GenomeID,Genome_ClusterID),by="Genome_ClusterID")%>%dplyr::rename(NRgenome_GenomeID=GenomeID.y)%>%dplyr::rename(GenomeID=GenomeID.x) # 9772, add NRgenome_GenomeID
CRBD_genome_meta_up2<-inner_join(CRBD_genome_meta_up,CRBD_genome_meta_up%>%filter(Rep_NR_Species=="Yes")%>%select(GenomeID,Species_ClusterID),by="Species_ClusterID")%>%dplyr::rename(RepSpecies_GenomeID=GenomeID.y)%>%dplyr::rename(GenomeID=GenomeID.x) #9772, add RepSpecies_GenomeID
CRBD_NRgenome<-CRBD_genome_meta_up2%>%filter(Rep_NR_Genomes=="Yes") #7531 CRBD NRgenome
CRBD_HQ_NRgenome<-CRBD_NRgenome%>%filter(Quality_level=="High Quality") #6109 CRBD HQ NRgenome
rownames(CRBD_HQ_NRgenome)<-CRBD_HQ_NRgenome$NRgenome_GenomeID
CRBD_RepSpecies<-CRBD_genome_meta_up2%>%filter(Rep_NR_Species=="Yes") # 3044 CRBD RepSpecies
CRBD_HQ_RepSpecies<-CRBD_RepSpecies%>%filter(Quality_level=="High Quality") # 2659 CRBD HQ RepSpecies
row.names(CRBD_HQ_RepSpecies)<-CRBD_HQ_RepSpecies$RepSpecies_GenomeID
CRBC_genomeID_list<-as.character((CRBD_HQ_NRgenome%>%filter(Source=="CRBC")%>%select(GenomeID))$GenomeID)

CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome%>%select(c(NRgenome_GenomeID,Phylum_gtdb,Class_gtdb,Order_gtdb,Family_gtdb,Genus_gtdb,Species_gtdb))
colnames(CRBD_HQ_NRgenome_taxonomy)<-c("NRgenome_GenomeID","Phylum","Class","Order","Family","Genus","Species")
CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome_taxonomy%>%mutate(PhyClass=if_else(Phylum=="p__Proteobacteria",Class,Phylum))
CRBD_HQ_NRgenome_taxonomy$PhyClass<-gsub("p__","",CRBD_HQ_NRgenome_taxonomy$PhyClass)
CRBD_HQ_NRgenome_taxonomy$PhyClass<-gsub("c__","",CRBD_HQ_NRgenome_taxonomy$PhyClass)
CRBD_HQ_NRgenome_taxonomy$Genus<-gsub("g__","",CRBD_HQ_NRgenome_taxonomy$Genus)
CRBD_HQ_NRgenome_taxonomy_PhyClass_size<-data.frame(CRBD_HQ_NRgenome_taxonomy%>%group_by(PhyClass)%>%mutate(PhyClass_size=n())%>%select(c(PhyClass,PhyClass_size))%>%unique()%>%arrange(desc(PhyClass_size))) # the order is not consistent with Ray's definition, but in order to keep the figure legend consistent, will define the dominant PhyClass just the same as Ray's
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/R_output/CRBD_HQ_NRgenome_PhyClass_size.pdf",height = 4,width = 6)

CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass,levels =unique(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass))
plot1 <- ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="HQ NRgenome number") +
  coord_cartesian(ylim = c(0, 100))
# Create the second segment of the plot
plot2 <- ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+labs(y="HQ NRgenome number") +
  coord_cartesian(ylim = c(1000, 2500)) 
# Combine the two segments into a single plot using facet_wrap
plot2+plot1+plot_layout(ncol=1)
#ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="HQ NRgenome number (log10)")
dev.off()

PhyClass_dominant=c("Alphaproteobacteria","Gammaproteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Myxococcota","Patescibacteria","Spirochaetota","Fibrobacterota","Chloroflexota","Verrucomicrobiota","Acidobacteriota") # Patescibacteria does not have HQ NRgenome
CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome_taxonomy%>%mutate(PhyClass_collapse=if_else(PhyClass%in%PhyClass_dominant,PhyClass,"Others"))
PhyClass_order<-c(PhyClass_dominant,"Others")
write.table(CRBD_HQ_NRgenome_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/CRBD_HQ_NRgenome_taxanomy_gtdb.txt",sep = "\t",quote=FALSE,row.names = FALSE)
rownames(CRBD_HQ_NRgenome_taxonomy)<-CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID # 6109 HQ RepSpecies
CRBD_HQ_NRgenome_taxonomy_PhyClass_size<-data.frame(CRBD_HQ_NRgenome_taxonomy%>%group_by(PhyClass_collapse)%>%summarise(PhyClass_size=n()))%>%arrange(desc(PhyClass_size))
HQ_NRgenome_taxonomy_up_sort<-inner_join(CRBD_HQ_NRgenome_taxonomy,CRBD_HQ_NRgenome_taxonomy_PhyClass_size)%>%arrange(desc(PhyClass_size))
rownames(HQ_NRgenome_taxonomy_up_sort)<-HQ_NRgenome_taxonomy_up_sort$NRgenome_GenomeID
HQ_NRgenome_taxonomy_up_sort$PhyClass_collapse<-factor(HQ_NRgenome_taxonomy_up_sort$PhyClass_collapse,levels = PhyClass_order[PhyClass_order!="Patescibacteria"]) #6109
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/R_output/CRBD_HQ_NRgenome_PhyClass_collapse_size.pdf",width=6,height = 4)
CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse,levels =unique(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse))
#ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="HQ NRgenome number")
plot1 <- ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="HQ NRgenome number") +
  coord_cartesian(ylim = c(0, 500))
# Create the second segment of the plot
plot2 <- ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+labs(y="HQ NRgenome number") +
  coord_cartesian(ylim = c(1000, 2500)) 
# Combine the two segments into a single plot using facet_wrap
plot2+plot1+plot_layout(ncol=1)
dev.off()

CRBD_RepSpecies_taxonomy_PhyClass_size<-data.frame(CRBD_RepSpecies%>%mutate(PhyClass=if_else(Phylum_gtdb=="p__Proteobacteria",Class_gtdb,Phylum_gtdb))%>%group_by(PhyClass)%>%mutate(PhyClass_size=n())%>%select(c(PhyClass,PhyClass_size))%>%unique()%>%arrange(desc(PhyClass_size)))
CRBD_RepSpecies_taxonomy_PhyClass_size$PhyClass<-factor(CRBD_RepSpecies_taxonomy_PhyClass_size$PhyClass,levels = unique(CRBD_RepSpecies_taxonomy_PhyClass_size$PhyClass))
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/R_output/CRBD_RepSpecies_PhyClass_collapse_size.pdf",width = 6,height = 4)
plot1 <- ggplot(CRBD_RepSpecies_taxonomy_PhyClass_size,aes(x=PhyClass,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="RepSpecies number") +
  coord_cartesian(ylim = c(0, 250))
# Create the second segment of the plot
plot2 <- ggplot(CRBD_RepSpecies_taxonomy_PhyClass_size,aes(x=PhyClass,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_blank(),axis.ticks.x =element_blank(),axis.title.x = element_blank())+labs(y="RepSpecies number") +
  coord_cartesian(ylim = c(700, 900))
# Combine the two segments into a single plot using facet_wrap
plot2+plot1+plot_layout(ncol=1)
#ggplot(CRBD_RepSpecies_taxonomy_PhyClass_size,aes(x=PhyClass,y=PhyClass_size))+geom_bar(stat = 'identity',width = 0.7,fill="#5a2341")+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(y="RepSpecies number")
dev.off()

PhyClass_col<-list(Color=taxa_color[c(1:12,15),]$Color)
names(PhyClass_col$Color)<-taxa_color[c(1:12,15),]$PhyClass

# BGC basic stats 

### ------- BGC_5Kbplus_per_NRgenome ------

### ------- BGC infor 
CRBD_NRgenome_BGC_5kb<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Sep_8_update/input_files/CRBD_NRgenome_5Kplus_BGC_cluster_info_up.txt",sep = "\t",header=FALSE,quote = "",col.names = c("GenomeID","BGC","ContigID","Start","End")) # 55285
length(unique(CRBD_NRgenome_BGC_5kb$GenomeID)) # 7278 genomes detected BGC out of 7531 genomes
length(unique(CRBD_NRgenome_BGC_5kb$BGC)) # 55285 BGCs 
CRBD_NRgenome_BGC_5kb<-CRBD_NRgenome_BGC_5kb%>%mutate(BGC_length=(End-Start+1))
summary(CRBD_NRgenome_BGC_5kb$BGC_length) # Min. 1stQu. Median Mean 3rdQu. Max., 5002   13506   21234   25971   35033  241352 
unique(CRBD_NRgenome_BGC_5kb$GenomeID[!CRBD_NRgenome_BGC_5kb$GenomeID%in%CRBD_NRgenome$NRgenome_GenomeID])
CRBD_NRgenome_BGC_5kb<-CRBD_NRgenome_BGC_5kb%>%dplyr::filter(!str_detect(GenomeID,"Iso_Wt_W139-2|Iso_Wt_WTC65-1|Iso_Wt_WTC65-2")) #55261
unique(CRBD_NRgenome_BGC_5kb$GenomeID[!CRBD_NRgenome_BGC_5kb$GenomeID%in%CRBD_NRgenome$NRgenome_GenomeID]) # 7275 NRgenome and 55261 BGCs
CRBD_HQ_NRgenome_BGC_5kb<-CRBD_NRgenome_BGC_5kb%>%filter(GenomeID%in%CRBD_HQ_NRgenome$NRgenome_GenomeID) #6087 HQ NRgenome out of 6109 HQ NRgenome has BGC potential, 48644

### ------- BigScape Class
CRBD_NRgenome_BGC_bigscape_Class<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Sep_8_update/input_files/CRBD_NRgenome_5Kplus_BGC_bigscapeClass_up.txt",sep="\t",col.names = c("BGC_ori","BigScape_Class")) # 55285 BGCs
CRBD_NRgenome_BGC_bigscape_Class$BGC<-paste(CRBD_NRgenome_BGC_bigscape_Class$BGC_ori,".gbk",sep="")
CRBD_NRgenome_BGC_bigscape_Class<-CRBD_NRgenome_BGC_bigscape_Class%>%select(!BGC_ori)
CRBD_NRgenome_BGC_bigscape_Class<-CRBD_NRgenome_BGC_bigscape_Class%>%dplyr::filter(!str_detect(BGC,"Iso_Wt_W139-2|Iso_Wt_WTC65-1|Iso_Wt_WTC65-2")) #55261

### -------merge BGC with Bigscape
CRBD_HQ_NRgenome_BGC_5kb_update_Class<-inner_join(CRBD_HQ_NRgenome_BGC_5kb,CRBD_NRgenome_BGC_bigscape_Class) # 48644 BGC across 6087 NRgenome 

### ------- BGC completeness
CRBD_BGC_completeness<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Sep_8_update/input_files/BGC_completeness.txt",sep = "\t",header = TRUE,quote = "")
summary(CRBD_BGC_completeness$BGC%in%CRBD_HQ_NRgenome_BGC_5kb_update_Class$BGC) # 48644 

### ------- BGC novelty
CRBD_BGC_novelty<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Sep_8_update/input_files/BGC_known_up.list",sep="\t",header = FALSE,col.names = c("GenomeID","BGC","Known_list"),quote = "")
CRBD_BGC_novelty$Known_list[CRBD_BGC_novelty$Known_list==""]="Unknown"
CRBD_BGC_novelty<-CRBD_BGC_novelty%>%mutate(Novelty=if_else(Known_list=="Unknown","Novel","Known"))%>%select(!c(GenomeID,Known_list))
summary(CRBD_BGC_novelty$BGC%in%CRBD_HQ_NRgenome_BGC_5kb_update_Class$BGC) # 48644, the genomeID of Iso_Wt_WTC89-2 is not changed to Iso_Wt_WTC89_2

### ------- BGC bigscape class + novelty + completeness
CRBD_HQ_NRgenome_BGC_5kb_update_Class_up<-data.frame(inner_join(inner_join(CRBD_HQ_NRgenome_BGC_5kb_update_Class,CRBD_BGC_completeness,by="BGC"),CRBD_BGC_novelty,by="BGC"))%>%dplyr::rename(Completeness=Complete) # 48644 BGC across 6087 HQ NRgenomes
CRBD_HQ_NRgenome_BGC_5kb_update_Class_up$BigScape_Class<-factor(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up$BigScape_Class,levels = BGC_class)


### ------  Bacic stats ----------

##### 1. 6087 genomes detected BGC out of 6109 genomes
##### 2. 55261 BGCs detected across 6087 genomes
##### 3. length summary: Min. 1stQu. Median Mean 3rdQu. Max., 5002   16236   21916   27431   40840  241352 
##### 4. Complete ration: 36569/48644=75.17679%
##### 5. Novel ratio: 27774/48644=57.09646%

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_NRgenome_BGC_5kbplus_completeness_vs_BGC_length.pdf",width = 6,height = 3)
ggplot(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up,aes(x=BGC_length/1000,color=Completeness,fill=Completeness))+geom_histogram(bins = 80,color="#787272")+theme_bw()+theme_bw()+labs(x="BGC length (Kb)",y="Frequency")+scale_fill_manual(values = c("#db6618","#a19e9c"))
dev.off()
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_NRgenome_BGC_5kbplus_completeness_vs_BGC_length_within_100kb.pdf",width = 6,height = 3)
ggplot(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up,aes(x=BGC_length/1000,color=Completeness,fill=Completeness))+geom_histogram(bins = 80,color="#787272")+xlim(5,100)+theme_bw()+labs(x="BGC length (Kb)",y="Frequency")+scale_fill_manual(values = c("#db6618","#a19e9c"))
dev.off()

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_NRgenome_BGC_5kbplus_BGC_length_between_BigScapeClasses.pdf",width = 6,height = 3)
ggplot(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up,aes(x=BGC_length/1000,color=Completeness,fill=BigScape_Class))+geom_histogram(bins = 80,color="#e1e3e3")+theme_bw()+theme_bw()+labs(x="BGC length (Kb)",y="Frequency")+scale_fill_manual(values = BGC_color)
dev.off()
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_NRgenome_BGC_5kbplus_BGC_length_within_100kb_between_BigScapeClass.pdf",width = 6,height = 3)
ggplot(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up,aes(x=BGC_length/1000,color=Completeness,fill=BigScape_Class))+geom_histogram(bins = 80,color="#e1e3e3",size=0.1)+xlim(5,100)+theme_bw()+labs(x="BGC length (Kb)",y="Frequency")+scale_fill_manual(values = BGC_color)
dev.off()

### ------- BGC frequency distribution ---
### ------- BGC frequency -----
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome<-data.frame(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up%>%select(c(GenomeID,BGC))%>%group_by(GenomeID)%>%mutate(BGC_Freq_perGenome=n())%>%select(!BGC)%>%unique())
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta<-full_join(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta$BGC_Freq_perGenome[is.na(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta$BGC_Freq_perGenome)]=0 # 6109 CRBD HQ NRgenome and 6087 has BGC and 22 does not has BGC

CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta<-CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta%>%arrange(PhyClass_collapse,Class,Order,Family,Genus,Species)
sum((CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta%>%filter(GenomeID%in%CRBC_genomeID_list))$BGC_Freq_perGenome) # 27893
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta$PhyClass_collapse,levels = as.character(PhyClass_order))
NRgenomeID_order<-as.character(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta$GenomeID)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_per_Genome_barplot.pdf",height = 4,width = 15)
ggplot(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta,aes(x=GenomeID,y=BGC_Freq_perGenome))+geom_bar(stat = "identity",width = 0.5,alpha=0.8,fill="#467d82")+facet_grid(~PhyClass_collapse,scales = "free", space = "free")+main_theme+theme(axis.text.x = element_blank(),strip.text = element_text(size=6))
dev.off()
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_per_Genome_fix_width_barplot.pdf",height = 4,width = 15)
ggplot(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta,aes(x=GenomeID,y=BGC_Freq_perGenome))+geom_bar(stat = "identity",width = 0.5,alpha=0.8,fill="#467d82")+facet_grid(~PhyClass_collapse,scales = "free", space = "fixed")+main_theme+theme(axis.text.x = element_blank(),strip.text = element_text(size=6))
dev.off()

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_per_Genome_boxplot_between_PhyClass.pdf",height =8,width = 8)
plotA<-ggplot(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta,aes(x=PhyClass_collapse,y=BGC_Freq_perGenome))+geom_boxplot(fill="#81b5b8",color='#377c80',outlier.alpha = 0.1,size=0.3)+main_theme+theme(axis.text.x = element_blank(),axis.title.x = element_blank())+labs(x="PhyClass",y="BGC richness")
CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_collapse,levels = PhyClass_order)
plotB<-ggplot(CRBD_HQ_NRgenome_taxonomy_PhyClass_size,aes(x=PhyClass_collapse,y=PhyClass_size))+geom_text(aes(label=PhyClass_size), position=position_dodge(width=0.9), vjust=-0.25)+geom_bar(stat = "identity",fill="#81b5b8")+main_theme+theme(axis.text.x = element_text(size=15,hjust = 1))+labs(x="PhyClass",y="NRgenome number")+ylim(0,2700)
plotA/plotB
dev.off()

### ------- BGC Freq vs genome size ----------

CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up<-inner_join(data.frame(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta%>%select(c(GenomeID,BGC_Freq_perGenome,PhyClass_collapse))%>%unique()),CRBD_genome_meta_up2 %>%select(c(GenomeID,Genome_size)))
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up$PhyClass_collapse,levels = PhyClass_order)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_vs_Genome_size_scatter_plot.pdf",width = 8,height = 6)
ggplot(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up,aes(x=Genome_size/1000000,y=BGC_Freq_perGenome,color=PhyClass_collapse,fill=PhyClass_collapse))+geom_point(size=0.5,alpha=0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+scale_color_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+labs(x="Genome size (Mb)",y="BGC frequence")
dev.off()

### ------- BGC Freq vs genome completeness ----------

CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up2<-inner_join(data.frame(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta%>%select(c(GenomeID,BGC_Freq_perGenome,PhyClass_collapse))%>%unique()),CRBD_genome_meta_up2 %>%select(c(GenomeID,Completeness)))
CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up2$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up2$PhyClass_collapse,levels = PhyClass_order)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_vs_Genome_completeness_scatter_plot.pdf",width = 8,height = 5)
ggplot(CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up2%>%filter(PhyClass_collapse!="Others"),aes(x=Completeness,y=BGC_Freq_perGenome,color=PhyClass_collapse,fill=PhyClass_collapse))+facet_wrap(vars(PhyClass_collapse),nrow = 4,ncol = 4)+geom_point(size=0.5)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+scale_color_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+labs(x="Genome completeness (%)",y="BGC frequence")
dev.off()

### ------- BGC encoding gene percentage -------

BGC_encoding_ratio_perGenome<-full_join(data.frame(CRBD_NRgenome_BGC_5kb%>%select(c(GenomeID,BGC,BGC_length))%>%unique()),CRBD_genome_meta_up2%>%filter(GenomeID%in%CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID )%>%select(c(GenomeID,Genome_size)))
BGC_encoding_ratio_perGenome$BGC_length[is.na(BGC_encoding_ratio_perGenome$BGC_length)]=0
BGC_encoding_ratio_perGenome<-data.frame(BGC_encoding_ratio_perGenome%>%group_by(GenomeID)%>%mutate(BGC_length_Sum=sum(BGC_length))%>%select(!c(BGC_length,BGC))%>%unique()%>%mutate(BGC_encoding_ratio_perGenome=BGC_length_Sum/Genome_size))
BGC_encoding_ratio_perGenome_taxa<-full_join(BGC_encoding_ratio_perGenome,CRBD_HQ_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,PhyClass_collapse)),by=c("GenomeID"="NRgenome_GenomeID"))
BGC_encoding_ratio_perGenome_taxa$PhyClass_collapse<-factor(BGC_encoding_ratio_perGenome_taxa$PhyClass_collapse,levels = PhyClass_order)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_encoding_percentage_vs_Genome_size_scatter_plot_individual_phylum.pdf",width = 10,height = 6)
ggplot(BGC_encoding_ratio_perGenome_taxa%>%filter(PhyClass_collapse!="Others"),aes(x=Genome_size/1000000,y=BGC_encoding_ratio_perGenome*100,color=PhyClass_collapse,fill=PhyClass_collapse))+geom_point(size=0.5)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+scale_color_manual(values = taxa_color[PhyClass_order,]$Color)+facet_wrap(vars(PhyClass_collapse),nrow = 4,ncol = 4)+theme_bw()+theme(axis.text.x = element_text(size=10))+labs(x="Genome size (Mb)",y="BGC encoding percentage (%)")
dev.off()


pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_encoding_percentage_vs_Genome_size_scatter_plot_combine_all_phyla.pdf",width = 7,height = 5)
ggplot(BGC_encoding_ratio_perGenome_taxa%>%filter(PhyClass_collapse!="Others"),aes(x=Genome_size/1000000,y=BGC_encoding_ratio_perGenome*100,color=PhyClass_collapse,fill=PhyClass_collapse))+geom_point(size=0.5)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+scale_color_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(size=10))+labs(x="Genome size (Mb)",y="BGC encoding percentage (%)")
dev.off()

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_genome_size_density_plot.pdf",height = 5,width = 5)
ggplot(data=BGC_encoding_ratio_perGenome_taxa%>%filter(PhyClass_collapse!="Others"),aes(x=Genome_size/1000000))+geom_density()+theme_bw()+labs(x="Genome size")
dev.off()

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_HQ_genome_BGC_encoding_precentage_density_plot.pdf",height = 5,width = 5)
ggplot(data=BGC_encoding_ratio_perGenome_taxa%>%filter(PhyClass_collapse!="Others"),aes(x=BGC_encoding_ratio_perGenome*100))+geom_density()+theme_bw()+labs(x="BGC encoding percentage")
dev.off()


### ------- BGC present ratio -------

CRBD_BGC_presence_stat<-CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome_meta_up2%>%mutate(Presence=if_else(BGC_Freq_perGenome>=1,"Presence","Absence"))%>%select(c(GenomeID,Presence))
CRBD_BGC_presence_stat_taxa<-data.frame(full_join(CRBD_BGC_presence_stat,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))%>%group_by(PhyClass_collapse,Presence)%>%mutate(Presence_perPhyClass=n())%>%select(c(PhyClass_collapse,Presence_perPhyClass,Presence))%>%unique())
sum(CRBD_HQ_NRgenome_taxonomy_PhyClass_size$PhyClass_size) # 6109
CRBD_BGC_presence_stat_taxa_up<-full_join(CRBD_BGC_presence_stat_taxa,CRBD_HQ_NRgenome_taxonomy_PhyClass_size)
CRBD_BGC_presence_stat_taxa_up$BGC_PresenceRatio_perPhyClass<-round(CRBD_BGC_presence_stat_taxa_up$Presence_perPhyClass/CRBD_BGC_presence_stat_taxa_up$PhyClass_size*100,digits = 2)
CRBD_BGC_presence_stat_taxa_up$PhyClass_collapse<-factor(CRBD_BGC_presence_stat_taxa_up$PhyClass_collapse,levels = PhyClass_order)


pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_BGC_presence_Genome_number_vs_BGC_richness.pdf",height = 8,width = 8)
plotC<-ggplot(CRBD_BGC_presence_stat_taxa_up,aes(x=PhyClass_collapse,y=Presence_perPhyClass,fill=Presence))+geom_bar(stat="identity",width = 0.8,alpha=0.9)+scale_fill_manual(values = c("#edaa32","#68967e"))+main_theme+theme(axis.text.x = element_text(size=15,hjust = 1))+labs(x="PhyClass",y="Genome number")+ylim(0,2700)
plotA/plotC
dev.off()


### ------- Novel BGC ratio ------

CRBD_NRgenome_BGC_5kb_NovelSum<-CRBD_HQ_NRgenome_BGC_5kb_update_Class_up%>%group_by(GenomeID,Novelty)%>%mutate(BGC_noveltySum=n())%>%select(c(GenomeID,Novelty,BGC_noveltySum))%>%unique()%>%filter(Novelty=="Novel") # 6927 out of 7275 BGC synthesis NRgenomes
CRBD_NRgenome_BGC_5kb_NovelSum<-full_join(data.frame(GenomeID=CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID),CRBD_NRgenome_BGC_5kb_NovelSum)
CRBD_NRgenome_BGC_5kb_NovelSum$Novelty[is.na(CRBD_NRgenome_BGC_5kb_NovelSum$Novelty)]="Novelty"
CRBD_NRgenome_BGC_5kb_NovelSum$BGC_noveltySum[is.na(CRBD_NRgenome_BGC_5kb_NovelSum$BGC_noveltySum)]=0

CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq<-full_join(CRBD_NRgenome_BGC_5kb_NovelSum%>%select(c(GenomeID,BGC_noveltySum)),CRBD_HQ_NRgenome_5bk_BGC_Freq_perGenome)
CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq$BGC_Freq_perGenome[is.na(CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq$BGC_Freq_perGenome)]=0
CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq<-CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq%>%mutate(NovelRatio=BGC_noveltySum/BGC_Freq_perGenome*100)
CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq_taxa<-full_join(CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq_taxa$PhyClass_collapse<-factor(CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq_taxa$PhyClass_collapse,levels = PhyClass_order)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_NRgenome_BGC_Freq_vs_Novelty_scatter_plot.pdf",width = 8,height = 5)
ggplot(CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq_taxa%>%filter(PhyClass_collapse!="Others"),aes(x=BGC_Freq_perGenome,y=BGC_noveltySum,color=PhyClass_collapse,fill=PhyClass_collapse))+facet_wrap(vars(PhyClass_collapse),nrow = 4,ncol = 4)+geom_point(size=0.5)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+scale_color_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+labs(x="Total BGC",y="Novel BGC")+xlim(0,60)+ylim(0,60)+geom_abline(intercept = 0,slope = 1,color='grey',lty='dashed')
dev.off()

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/basic_stats/CRBD_BGC_PresenceNum_vs_BGC_richness_vs_Novelty_Ratio.pdf",width = 6,height = 10)
plotD<-ggplot(CRBD_NRgenome_BGC_5kb_NovelSum_BGC_Freq_taxa,aes(x=PhyClass_collapse,y=NovelRatio))+geom_boxplot(fill="#81b5b8",color='#377c80',outlier.alpha = 0.1,size=0.3)+scale_fill_manual(values = c("#edaa32","#68967e"))+main_theme+theme(axis.text.x = element_blank())+labs(x="PhyClass",y="Novelty")
plotA/plotD/plotC
dev.off()

# BGC composition - stacked barplot

CRBD_HQ_NRgenome_BGC_5kb_update_Class_up2<-CRBD_HQ_NRgenome_BGC_5kb_update_Class_up%>%select(c(GenomeID,BGC,BigScape_Class,Novelty))
CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa<-inner_join(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up2,CRBD_HQ_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,PhyClass_collapse))%>%dplyr::rename(GenomeID=NRgenome_GenomeID))
CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa_Rwrite<-inner_join(CRBD_HQ_NRgenome_BGC_5kb_update_Class_up2,CRBD_HQ_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,PhyClass_collapse,Class,Order,Family,Genus))%>%dplyr::rename(GenomeID=NRgenome_GenomeID))
write.table(CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa_Rwrite,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/R_write_CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa.txt",sep="\t",quote=FALSE,row.names = FALSE)
CRBD_NRgenome_BGCSum_perPhyClass<-CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa%>%group_by(PhyClass_collapse)%>%mutate(BGCsum_perPhyClass=n())%>%select(c(PhyClass_collapse,BGCsum_perPhyClass))%>%unique()
CRBD_NRgenome_BGCSum_perBigscapeClass_perPhyclass<-CRBD_NRgenome_BGC_5kb_update_bigscapeClass_taxa%>%group_by(PhyClass_collapse,BigScape_Class)%>%mutate(BGCSum_perBigscapeClass_perPhyClass=n())%>%select(c(PhyClass_collapse,BigScape_Class,BGCSum_perBigscapeClass_perPhyClass))%>%unique()
CRBD_NRgenome_BGCSum_perBigscapeClass_RA<-inner_join(CRBD_NRgenome_BGCSum_perPhyClass,CRBD_NRgenome_BGCSum_perBigscapeClass_perPhyclass)%>%mutate(BigscapeClass_RA=BGCSum_perBigscapeClass_perPhyClass/BGCsum_perPhyClass*100)
CRBD_NRgenome_BGCSum_perBigscapeClass_RA$PhyClass_collapse<-factor(CRBD_NRgenome_BGCSum_perBigscapeClass_RA$PhyClass_collapse,levels=PhyClass_order)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/BGC_composition_stacked_barplot.pdf",height = 4,width = 5)
ggplot(CRBD_NRgenome_BGCSum_perBigscapeClass_RA,aes(x=PhyClass_collapse,y=BigscapeClass_RA,fill=BigScape_Class))+geom_bar(stat='identity',colour="NA",size=0.1)+labs(x="PhyClass",y="Relative Abundance \n")+scale_fill_manual(values=BGC_color)+guides(fill=guide_legend(ncol=1))+main_theme+theme(legend.text = element_text(size=6,face = 'plain'),legend.title = element_text(size=8,face = 'plain'),legend.position = "right",axis.title.x= element_blank(),axis.title.y= element_text(size=8),axis.text.y = element_text(face = 'plain',size=6),axis.text.x = element_text(face = 'plain',size=8),strip.text = element_text(size=8,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"))
dev.off()

# generate the mapping between Genus and RepSpecies Genome

## extract the best RepSpecies as the representative genome for each Genus

CRBD_HQ_RepSpecies_sub<-data.frame(CRBD_HQ_RepSpecies%>%select(RepSpecies_GenomeID,Genus_gtdb,Score)%>%group_by(Genus_gtdb)%>%mutate(MaxScore=max(Score))%>%filter(Score==MaxScore))  #442 genus with g__ unclassified
summary(CRBD_HQ_RepSpecies_sub$Genus_gtdb=="g__") # 441 genus and 1 g__ unclassified
CRBD_HQ_RepSpecies_sub$Genus_gtdb<-gsub("g__","",CRBD_HQ_RepSpecies_sub$Genus_gtdb)
CRBD_HQ_RepSpecies_sub<-CRBD_HQ_RepSpecies_sub%>%filter(Genus_gtdb!="")%>%dplyr::rename(Genus=Genus_gtdb)
Genus_RepSpecies_map<-CRBD_HQ_RepSpecies_sub%>%select(RepSpecies_GenomeID,Genus)

## subset RepSpecies tree to genus representative ones

## ---- IRBC rep species tree ---

library(ggtree)
library(ape)
CRBD_RepSpecies_tree<-read.tree("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PhyloTree/CRBD_3044.unrooted.tree")
tree_order<-fortify(CRBD_RepSpecies_tree)%>%filter(isTip=="TRUE")%>%arrange(y) # 3044 RepSpecies
RepSpGenome_order_intree<-as.character(tree_order$label) # 3044
RepSpecies_order<-RepSpGenome_order_intree
nodes_to_trim<-RepSpecies_order[!RepSpecies_order%in%Genus_RepSpecies_map$RepSpecies_GenomeID] # 2603
CRBD_HQ_RepGenus_tree<-drop.tip(CRBD_RepSpecies_tree, nodes_to_trim)
length((fortify(CRBD_HQ_RepGenus_tree)%>%filter(isTip=="TRUE")%>%arrange(y))$label) # 441 CRBD HQ RepGenus
HQ_RepGenus_tree_order<-as.character((fortify(CRBD_HQ_RepGenus_tree)%>%filter(isTip=="TRUE")%>%arrange(y))$label) # 441
write.tree(CRBD_HQ_RepGenus_tree,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PhyloTree/CRBD_HQ_441_RepGenus_RepSpecies_unrooted.tree")
write.tree(CRBD_HQ_RepGenus_tree,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/CRBD_HQ_441_RepGenus_RepSpecies_unrooted.tree")

# BGC composition RepGenus Mean - phylogenetic tree

## --- HQ RepGenus mean ----

CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass<-data.frame(CRBD_HQ_NRgenome_BGC_5kb_BigScapeClass_add_taxa%>%select(c(NRPS,Others,`PKS-NRP_Hybrids`,PKSI,PKSother,RiPPs,Saccharides,Terpene,Genus))%>%group_by(Genus)%>%mutate_at(vars(-Genus),function(x) sum(x)/sum(x>=0))%>%unique())
CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up<-inner_join(CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass,Genus_RepSpecies_map)
row.names(CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up)<-CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up$RepSpecies_GenomeID
CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up<-CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up[HQ_RepGenus_tree_order,]

## ----- Tree annot ------

write.csv(CRBD_HQ_RepGenus_BGC_5kb_BigScapeClass_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/CRBD_HQ_NRgenome_BGC_RepGenus_BigScape_Class_RA.csv",row.names = FALSE,quote = FALSE)
source("/Users/fangliu/Documents/IGDB_Bai_lab/Script_backup/table2itol/table2itol-master/table2itol.R")
#Apparently this script is running in interactive mode. You could now generate iTOL files by setting some 'infiles' variable to a vector of file names and then calling:create_itol_files(infiles)
setwd("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/")
create_itol_files(infiles = "CRBD_HQ_NRgenome_BGC_RepGenus_BigScape_Class_RA.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")

# CRBD vs MiBIG novel GCF - Novel GCF composition and PhyClass composition
# NOTE!!! since the distance file is very large (21G and subset to 893M), I decided to have this section done on meta server instead of my local compuster


# GCF richness per genome, based on GCF results

### ----- GCF membership and novelty
GCF_membership<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BiGSlice/CRBD_HQ_NRgenome_BGC_GCF_membership.txt",sep="\t",header = TRUE)
CRBD_HQ_NRgenome_vs_BigFam_novel<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BiGSlice/R_write_CRBD_vs_BigFam_mean_cos_dist_perGCF.txt",sep = "\t",header = TRUE)
BGC_name_vs_id<-read.csv("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BiGSlice/BGC_name_vs_bgc_id_map.txt",sep=",",header = TRUE)
BGC_name_vs_id<-BGC_name_vs_id%>%dplyr::rename(bgc_id=id)%>%dplyr::rename(BGC=orig_filename)

CRBD_HQ_NRgenome_BGC_GCF_vs_BigFam_cosdist<-inner_join(inner_join(BGC_name_vs_id,GCF_membership),CRBD_HQ_NRgenome_vs_BigFam_novel%>%select(c(gcf_id,min_cosdist_GCF_mean)))
CRBD_HQ_NRgenome_BGC_GCF_vs_BigFam_cosdist_up<-inner_join(CRBD_HQ_NRgenome_BGC_GCF_vs_BigFam_cosdist,CRBD_HQ_NRgenome_BGC_5kb_update_Class_up%>%select(c(BGC,GenomeID)))%>%dplyr::rename(GCF=gcf_id) #48644 BGC clustered into 12865

## GCF sum perNRgenome

CRBD_HQ_NRgenome_GCF_sum_perNRgenome<-CRBD_HQ_NRgenome_BGC_GCF_vs_BigFam_cosdist_up%>%select(c(GenomeID,GCF))%>%unique()%>%group_by(GenomeID)%>%mutate(GCF_sum_perNRgenome=n())%>%select(c(GenomeID,GCF_sum_perNRgenome))%>%unique() # 6087 genomes
CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa<-data.frame(full_join(CRBD_HQ_NRgenome_GCF_sum_perNRgenome,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID")))
CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa$GCF_sum_perNRgenome[is.na(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa$GCF_sum_perNRgenome)]=0
CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa$PhyClass_collapse<-factor(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa$PhyClass_collapse,levels = PhyClass_order)
R_write_CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa<-CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa%>%arrange(PhyClass_collapse,Class,Order,Family,Genus,desc(GCF_sum_perNRgenome))
write.table(R_write_CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/R_write_CRBD_HQ_NRgenome_GCF_0.65_richnes_GCF_sum_perNRgenome.txt",sep = "\t",row.names = FALSE,quote = FALSE)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/CRBD_HQ_NRgenome_GCF_cos.2_richness_and_BGC_presence_or_absence.pdf",width = 8,height = 8)
plotE<-ggplot(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa,aes(x=PhyClass_collapse,y=GCF_sum_perNRgenome))+geom_boxplot(fill="#81b5b8",color='#377c80',outlier.alpha = 0.1,size=0.3)+main_theme+theme(axis.text.x = element_blank(),axis.title.x = element_blank())+labs(x="PhyClass",y="GCF richness")
plotE/plotC
dev.off()

Genus_size<-data.frame(CRBD_HQ_NRgenome_taxonomy%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%group_by(Genus)%>%mutate(Genus_size=n())%>%select(c(Genus,Genus_size))%>%unique()) #441
Genus_size<-inner_join(Genus_size,Genus_RepSpecies_map)
write.table(Genus_size,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/R_write_CRBD_HQ_NRgenome_RepGenus_size.txt",sep = "\t",row.names = FALSE,quote = FALSE)

## --- Actinobacteria GCF dissect on genus level -----

CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria<-CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa%>%filter(PhyClass_collapse=="Actinobacteriota")%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%arrange(desc(GCF_sum_perNRgenome))
CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Genus<-factor(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Genus,levels=unique(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Genus))
CRBD_HQ_NRgenome_Actinobacteria_genus_size<-Genus_size%>%filter(Genus%in%unique(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Genus))
CRBD_HQ_NRgenome_Actinobacteria_genus_size$Genus<-factor(CRBD_HQ_NRgenome_Actinobacteria_genus_size$Genus,levels =unique(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Genus))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/Actinobacteria_GCF_richness_perNRgenome_versus_Genus_size.pdf",width = 12,height = 8)
plot_Actinobacteria_GCF_richness<-ggplot(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria,aes(x=Genus,y=GCF_sum_perNRgenome))+geom_boxplot(fill="#81b5b8",color='#377c80',outlier.alpha = 0.1,size=0.3)+main_theme+theme(axis.text.x = element_blank(),axis.title.x = element_blank())+labs(x="Genus",y="GCF richness")
plot_Actinobacteria_Genus_size<-ggplot(CRBD_HQ_NRgenome_Actinobacteria_genus_size,aes(x=Genus,y=Genus_size))+geom_bar(fill="#81b5b8",color='#377c80',stat = 'identity',width = 0.7)+main_theme+theme(axis.text.x = element_text(family="Arial",size=10,angle = 90,hjust = 0,vjust = 0),axis.text.y = element_text(family="Arial",size=10))+labs(x="Genus",y="Genus size")
plot_Actinobacteria_GCF_richness/plot_Actinobacteria_Genus_size
dev.off()

## --- Actinobacteria GCF dissect on family level -----

CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF<-CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria%>%group_by(Family)%>%mutate(Family_meanGCF=mean(GCF_sum_perNRgenome))%>%select(c(Family,Family_meanGCF))%>%unique()%>%as.data.frame()%>%arrange(desc(Family_meanGCF)) # set up the order of family based on their GCF mean, will be used for box-dot plot
CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family<-factor(CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family,levels = as.character(unique(CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family)))
CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Family<-factor(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria$Family,levels =as.character(unique(CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family)))
Actinobacteria_family_meanGCF <-ggplot(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_Actinobacteria,aes(x=Family,y=GCF_sum_perNRgenome))+geom_boxplot(fill="#81b5b8",color='#377c80',outlier.alpha = 0.1,size=0.3)+main_theme+labs(x="Family",y="GCF richness")

Family_size_sub<-CRBD_HQ_NRgenome_taxonomy%>%filter(PhyClass_collapse=="Actinobacteriota")%>%group_by(Family)%>%mutate(Family_size=n())%>%select(c(Family,Family_size))%>%unique()%>%as.data.frame()%>%filter(Family%in%CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family)
Family_size_sub$Family<-factor(Family_size_sub$Family,levels =as.character(unique(CRBD_HQ_NRgenome_Actinobacteria_family_meanGCF$Family)))
Actinobacteria_family_size<-ggplot(Family_size_sub,aes(x=Family,y=Family_size))+geom_bar(fill="#81b5b8",color='#377c80',stat = 'identity',width = 0.7)+main_theme+theme(axis.text.x = element_text(family="Arial",size=10,angle = 90,hjust = 1,vjust = 0),axis.text.y = element_text(family="Arial",size=10))+labs(x="Family",y="Family size")

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/Actinobacteria_family_meanGCF_perNRgenome_versus_Family_size.pdf",width = 12,height = 8)
Actinobacteria_family_meanGCF/Actinobacteria_family_size
dev.off()

## Alphaproteobacteria
## Gammaproteobacteria
## Bacteroidota
## Firmicutes
## Myxococcota


## --- GCF Mean perGenus ------

CRBD_HQ_NRgenome_GCF_richness_GenusMean<-data.frame(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%group_by(Genus)%>%mutate(GCF_Mean_perGenus=mean(GCF_sum_perNRgenome))%>%select(c(Genus,GCF_Mean_perGenus))%>%unique())
summary(is.na(CRBD_HQ_NRgenome_GCF_richness_GenusMean$GCF_Mean_perGenus)) # five of the Genus does not has BGC, Nocardioides, 
CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa<-left_join(CRBD_HQ_NRgenome_GCF_richness_GenusMean,CRBD_HQ_NRgenome_taxonomy%>%select(PhyClass_collapse,Class,Genus)%>%unique())%>%arrange(PhyClass_collapse,desc(GCF_Mean_perGenus)) #441
CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa_up<-left_join(Genus_RepSpecies_map,CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa)
R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa_up_sort<-CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa_up%>%arrange(PhyClass_collapse,desc(GCF_Mean_perGenus))
write.table(R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa_up_sort,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/R_write_CRBD_HQ_NRgenome_GCF_Mean_perGenus.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMean_taxa_up_sort,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/R_write_CRBD_HQ_NRgenome_GCF_Mean_perGenus.txt",sep = "\t",quote = FALSE,row.names = FALSE)

## --- GCF Median perGenus ------

CRBD_HQ_NRgenome_GCF_richness_GenusMedian<-data.frame(CRBD_HQ_NRgenome_GCF_sum_perNRgenome_taxa%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%group_by(Genus)%>%mutate(GCF_Median_perGenus=median(GCF_sum_perNRgenome))%>%select(c(Genus,GCF_Median_perGenus))%>%unique())
summary(is.na(CRBD_HQ_NRgenome_GCF_richness_GenusMedian$GCF_Median_perGenus)) # five of the Genus does not has BGC, Nocardioides, 
CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa<-left_join(CRBD_HQ_NRgenome_GCF_richness_GenusMedian,CRBD_HQ_NRgenome_taxonomy%>%select(PhyClass_collapse,Class,Genus)%>%unique())%>%arrange(PhyClass_collapse,desc(GCF_Median_perGenus)) #441
CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa_up<-left_join(Genus_RepSpecies_map,CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa)
R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa_up_sort<-CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa_up%>%arrange(PhyClass_collapse,desc(GCF_Median_perGenus))
write.table(R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa_up_sort,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/GCF_richness/R_write_CRBD_HQ_NRgenome_GCF_Median_perGenus.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(R_write_CRBD_HQ_NRgenome_GCF_richness_GenusMedian_taxa_up_sort,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/R_write_CRBD_HQ_NRgenome_GCF_Median_perGenus.txt",sep = "\t",quote = FALSE,row.names = FALSE)

# Novel GCF per Genus (CRBD vs BigFam)

CRBD_HQ_NRgenome_novel_GCF<-CRBD_HQ_NRgenome_BGC_GCF_vs_BigFam_cosdist_up%>%mutate(Novelty=if_else(min_cosdist_GCF_mean>0.2,"Novel","Known")) # novel GCF are defined by mean min cos distance between CRBD vs BigFam
CRBD_HQ_NRgenome_novel_GCF_taxa<-inner_join(CRBD_HQ_NRgenome_novel_GCF,CRBD_HQ_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,PhyClass_collapse,Genus)),by=c("GenomeID"="NRgenome_GenomeID")) #48643
CRBD_HQ_NRgenome_GCF_sum_perGenome<-data.frame(CRBD_HQ_NRgenome_novel_GCF_taxa%>%select(c(GCF,GenomeID,Genus))%>%unique()%>%group_by(GenomeID)%>%mutate(GCF_sum_perGenome=n())%>%select(GenomeID,GCF_sum_perGenome,Genus)%>%unique()) # 6087 
CRBD_HQ_NRgenome_novel_GCF_sum_perGenome<-data.frame(CRBD_HQ_NRgenome_novel_GCF_taxa%>%filter(Novelty=="Novel")%>%select(c(GCF,GenomeID,Genus))%>%unique()%>%group_by(GenomeID)%>%mutate(Novel_GCF_sum_perGenome=n())%>%select(GenomeID,Novel_GCF_sum_perGenome,Genus)%>%unique()) # 4243 only
CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome<-full_join(CRBD_HQ_NRgenome_GCF_sum_perGenome,CRBD_HQ_NRgenome_novel_GCF_sum_perGenome%>%select(c(GenomeID,Novel_GCF_sum_perGenome)))
CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome$Novel_GCF_sum_perGenome[is.na(CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome$Novel_GCF_sum_perGenome)]=0
CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome<-CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome%>%mutate(Novel_GCF_ratio_perGenome=Novel_GCF_sum_perGenome/GCF_sum_perGenome*100) # 6087
CRBD_HQ_NRgenome_GCF_novelty_perGenus<-data.frame(CRBD_HQ_NRgenome_novel_GCF_ratio_perGenome%>%group_by(Genus)%>%mutate(Novelty_perGenus=mean(Novel_GCF_ratio_perGenome))%>%select(Genus,Novelty_perGenus)%>%unique()) # 444 genus
R_write_CRBD_HQ_NRgenome_GCF_novelty_perGenus_up<-inner_join(CRBD_HQ_NRgenome_GCF_novelty_perGenus,Genus_RepSpecies_map)%>%select(c(RepSpecies_GenomeID,Genus,Novelty_perGenus)) # 441 genus matched with Tree
write.table(R_write_CRBD_HQ_NRgenome_GCF_novelty_perGenus_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/R_wrire_CRBD_HQ_NRgenome_novel_GCF_vs_BiG_FAM_ratio.txt",sep = "\t",row.names = FALSE,quote = FALSE)

## ------ Tree annot -------

source("/Users/fangliu/Documents/IGDB_Bai_lab/Script_backup/table2itol/table2itol-master/table2itol.R")
#Apparently this script is running in interactive mode. You could now generate iTOL files by setting some 'infiles' variable to a vector of file names and then calling:create_itol_files(infiles)
setwd("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/BGC/Nov_6_HQ_NRgenome/BGC_composition_and_Tree/Tree_annot/")
create_itol_files(infiles = "R_wrire_CRBD_HQ_NRgenome_novel_GCF_vs_BiG_FAM_ratio.txt",identifier = "RepSpecies_GenomeID",label ="RepSpecies_GenomeID",separator = "\t")
create_itol_files(infiles = "R_write_CRBD_HQ_NRgenome_RepGenus_size.txt",identifier = "RepSpecies_GenomeID",label ="RepSpecies_GenomeID",separator = "\t")
create_itol_files(infiles = "R_write_CRBD_HQ_NRgenome_GCF_Median_perGenus.txt",identifier = "RepSpecies_GenomeID",label ="RepSpecies_GenomeID",separator = "\t",double.to.bars = TRUE)
```

