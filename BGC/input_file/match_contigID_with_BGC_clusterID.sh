# ----- link BGC cluster file to current directory ----------                                                             
ln -s /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/cluster_info.txt ./
wc -l CRBD_NRgenomeID_list #7531 CRBD_NRgenomeID_list
grep -w -F -f CRBD_NRgenomeID_list cluster_info.txt | cut -f 1 | sort | uniq | wc -l # 7393 has BGC, use below script to check whether those genomes does not has BGC or not
grep -v -w -F -f <(cut -f 1 cluster_info.txt| sort | uniq) CRBD_NRgenomeID_list  > temp
for file in $(cat temp); do echo $file; ls /mnt/m2/dairui/project/binning/MAG_finalization/all/14_BGC/sep_genome/"$file"/; done # it shows those genomes indeed lack of region.gbk file (BGC file). However, there is one genome has inconsistent genomeID between Ray BGC folder and CRBD genomeID list except Iso_Wt_WTC89_2 

## NOTE!!! # Iso_Wt_WTC89_2 is the genomeID within CRBD genomeID list, however, Ray BGC folder name changed to Iso_Wt_WTC89-2, so I need to update the BGC file by changing Iso_Wt_WTC89-2 to Iso_Wt_WTC89_2 -------
sed 's,Iso_Wt_WTC89-2,Iso_Wt_WTC89_2,g' cluster_info.txt > cluster_info_up.txt # remember to change the code when soft link BGC for bigscape analysis. Another three "Iso_Wt_W139-2"  "Iso_Wt_WTC65-1" "Iso_Wt_WTC65-2" genomeID is not consistent between BGC and taxonomy, since they are not included in CRBD NRgenomeID, I did not change it.

## NOTE!!! 274 rows within the cluster_info_up.txt files lack of contigID, when I checked Ray BGC output, it seemds the DEFINITION row is empty for those rows 
cat <(grep -v -P '\t\t' cluster_info_up.txt)  <(paste -d '\t'  <(grep -P '\t\t' cluster_info_up.txt | cut -f 1,2) <(grep -P '\t\t' cluster_info_up.txt | cut -f 2 | sed 's,.region.*.gbk,,g') <(grep -P '\t\t' cluster_info_up.txt | cut -f 4-)) > cluster_info_up2.txt # 134489


# -----------------  below the cluster BGC gbk column should not be changed since later this gbk list will be used to soft link for bigscape ect. ----- In order to match contig length file with BGC file contigID, Since lots of Pub contigID within BGC cluster file is pruned to be shorter than those within contig length file, I dicided to first split BGC cluster file into two section: Iso, MAG and Pub


# ----- in order to merge contig file within BGC contig information, I first subset cluster_info_up.txt into Iso, MAG and Pub, three parts ------

grep '^Iso' cluster_info_up2.txt > Iso_cluser_info_up2.txt
grep '^MAG' cluster_info_up2.txt > MAG_cluster_info_up2.txt
grep '^Pub' cluster_info_up2.txt > Pub_cluster_info_up2.txt

### ---- check the contigID consistence of Isolate -------

cut -f 3 Iso_cluser_info_up2.txt | sort | uniq | wc -l # 27020
grep -w -F -f <(cut -f 3 Iso_cluser_info_up2.txt | sort | uniq) <(cut -f 2 all_14242_genomes_contig_length.txt) | sort | uniq | wc -l #27000 # 20 of the contigs were absent within all_14242_genomes_contig_length.txt.
grep -v -w -F -f <(cut -f 2 all_14242_genomes_contig_length.txt) <(cut -f 3 Iso_cluser_info_up2.txt | sort | uniq)
```
Iso_Wt_W139-2___NODE_11_length_202916_cov_45.299331
Iso_Wt_W139-2___NODE_1_length_905099_cov_52.652907
Iso_Wt_W139-2___NODE_2_length_790113_cov_55.031673
Iso_Wt_W139-2___NODE_3_length_769478_cov_46.183656
Iso_Wt_W139-2___NODE_4_length_712422_cov_49.601971
Iso_Wt_W139-2___NODE_5_length_554489_cov_47.653200
Iso_Wt_WTC65-1___NODE_11_length_163689_cov_115.193507
Iso_Wt_WTC65-1___NODE_15_length_129870_cov_127.853051
Iso_Wt_WTC65-1___NODE_21_length_54006_cov_109.983645
Iso_Wt_WTC65-1___NODE_2_length_760481_cov_104.418347
Iso_Wt_WTC65-1___NODE_5_length_324423_cov_113.261650
Iso_Wt_WTC65-1___NODE_7_length_285605_cov_109.193182
Iso_Wt_WTC65-1___NODE_9_length_217648_cov_110.008158
Iso_Wt_WTC65-2___NODE_12_length_163690_cov_91.594922
Iso_Wt_WTC65-2___NODE_15_length_145368_cov_86.490904
Iso_Wt_WTC65-2___NODE_17_length_141950_cov_86.133852
Iso_Wt_WTC65-2___NODE_18_length_129872_cov_101.168342
Iso_Wt_WTC65-2___NODE_1_length_760502_cov_81.976504
Iso_Wt_WTC65-2___NODE_26_length_65522_cov_85.383467
Iso_Wt_WTC65-2___NODE_5_length_324423_cov_90.037004
```
### ---- to correct the above Isolate genomeID within the contig file since the genomeID within cluster file is consistent with that of the CRBD NRgenomeID-------

sed 's,Iso_Wt_W139_2,Iso_Wt_W139-2,g' all_14242_genomes_contig_length.txt | sed 's,Iso_Wt_WTC65_1,Iso_Wt_WTC65-1,g' | sed 's,Iso_Wt_WTC65_2,Iso_Wt_WTC65-2,g' > all_14242_genomes_contig_length_up.txt # then after double check the contig consistency between contig file and Isolate cluster file
cut -f 3 Iso_cluser_info_up2.txt | sort | uniq | wc -l # 27020
grep -w -F -f <(cut -f 3 Iso_cluser_info_up2.txt | sort | uniq) <(cut -f 2 all_14242_genomes_contig_length_up.txt) | sort | uniq | wc -l # 27020

### ------ check the contigID consistence of MAG --------

cut -f 3 MAG_cluster_info_up2.txt | sort | uniq | wc -l # 28838
grep -w -F -f <(cut -f 3 MAG_cluster_info_up2.txt | sort | uniq) <(cut -f 2 all_14242_genomes_contig_length.txt) | sort | uniq | wc -l #27781
grep -v -w -F -f <(cut -f 2 all_14242_genomes_contig_length.txt) <(cut -f 3 MAG_cluster_info_up2.txt | sort | uniq)

```
MAG_Mz_BJSZSb019___k127_612132
MAG_Mz_BJSZSb020___k127_105112
MAG_Mz_BJSZSb020___k127_110550
MAG_Mz_BJszSb020___k127_115916
MAG_Mz_BJszSb020___k127_153992
MAG_Mz_BJszSb020___k127_181579
MAG_Mz_BJszSb020___k127_229301
MAG_Mz_BJSZSb020___k127_238173
MAG_Mz_BJszSb020___k127_259114
MAG_Mz_BJszSb020___k127_309593
...
```
# ----- correct the MAG cluster contig column only since the NRgenomeID is correct and used BJszY21 and BJszY22 ----

paste -d '\t' <(cut -f 1,2 MAG_cluster_info_up2.txt) <(cut -f 3 MAG_cluster_info_up2.txt | sed 's,MAG_Mz_BJsz,MAG_Mz_BJszY21,g' | sed 's,MAG_Mz_BJSZ,MAG_Mz_BJszY22,g') <(cut -f 4- MAG_cluster_info_up2.txt) > MAG_cluster_info_up3.txt 
cut -f 3 MAG_cluster_info_up3.txt | sort | uniq | wc -l # 28838
grep -w -F -f <(cut -f 3 MAG_cluster_info_up3.txt | sort | uniq) <(cut -f 2 all_14242_genomes_contig_length_up.txt) | sort | uniq | wc -l # now they are just the same

# ----- solve the Pub contigID inconsistency between contig length file and BGC cluster file --------

wc -l  Pub_cluster_info_up2.txt # 60987 
paste -d '\t' <(cut -f 1,2 Pub_cluster_info_up2.txt) <(cut -f 3 Pub_cluster_info_up2.txt | sed 's,__,___,g' | sed 's,____,___,g') <(cut -f 4- Pub_cluster_info_up2.txt) > Pub_cluster_info_up3.txt  #60987 since the contigID use __ instead of ___, first correct this section, only change the contigID column and the GenomeID column is correct and the gbk column should not be changed since later analysis will need this list

### ------- whole word match contig list within cluste file --------
grep -w -F -f <(cut -f 3 Pub_cluster_info_up3.txt | sort | uniq) all_14242_genomes_contig_length.txt | cut -f 2 > Pub_cluster_word_match_list_in_contig_length # 10091
grep -w -F -f Pub_cluster_word_match_list <(cut -f 3 Pub_cluster_info_up3.txt | sort | uniq)  > Pub_cluster_word_match_list_exactly_in_cluster_file # 8932
grep -v -w -F -f Pub_cluster_word_match_list_exactly_in_cluster_file  Pub_cluster_word_match_list_in_contig_length | sort | uniq # it turned out some of the word match using cluster contigID as the query fished some of the contigID that are not full cover of the contigID within contig length file. Since some of the word match actualy followed by double round bracket 
```
B022DRAFT_2524025317.11OxalobacteraceaebacteriumJGI0001004-K23(contaminationscreened):B022DRAFT_2524025317.11
B022DRAFT_2524025322.16OxalobacteraceaebacteriumJGI0001004-K23(contaminationscreened):B022DRAFT_2524025322.16
B023DRAFT_2524025278.6OxalobacteraceaebacteriumJGI0001004-J12(contaminationscreened):B023DRAFT_2524025278.6
B136DRAFT_2524026109.118Rhizobiumsp.JGI0001005-K05(contaminationscreened):B136DRAFT_2524026109.118
...
...
```
grep 'B022DRAFT_2524025317.11OxalobacteraceaebacteriumJGI0001004-K23' Pub_cluster_info_up3.txt 
```
Pub_2528768082	B022DRAFT_2524025317.11.region001.gbk	B022DRAFT_2524025317.11OxalobacteraceaebacteriumJGI0001004-K23	26964	48749	other	hserlactone
```
grep "B022DRAFT_2524025317.11OxalobacteraceaebacteriumJGI0001004-K23" all_14242_genomes_contig_length_up.txt 
```
Pub_2528768082	B022DRAFT_2524025317.11OxalobacteraceaebacteriumJGI0001004-K23(contaminationscreened):B022DRAFT_2524025317.11	63840
```
cp Pub_cluster_word_match_list_exactly_in_cluster_file Pub_cluster_exactly_word_match_list # 8932
grep -v -w -F -f Pub_cluster_exactly_word_match_contig_list <(cut -f 3 Pub_cluster_info_up3.txt | sort | uniq) > Pub_cluster_partial_match_list # 31979
cut -f 3 Pub_cluster_info_up3.txt | sort | uniq | wc -l #40911

###  ----------- seperate contig length file based on exact word match and partial match -----------

cat Pub_cluster_exactly_word_match_list <(cut -f 3 Iso_cluser_info_up2.txt | sort | uniq) <(cut -f 3 MAG_cluster_info_up3.txt | sort | uniq) > Iso_MAG_Pub_exact_word_match_within_cluster_list #64790
grep -w -F -f Iso_MAG_Pub_exact_word_match_within_cluster_list all_14242_genomes_contig_length_up.txt > all_14242_genomes_contig_length_up_exact_word_match_subset.txt #64790 
wc -l  Pub_cluster_partial_match_list #31979
grep -F -f Pub_cluster_partial_match_list all_14242_genomes_contig_length_up.txt > all_14242_genomes_contig_length_up_partial_match_subset.txt # 31978
grep -v -F -f <(grep -o -F -f Pub_cluster_partial_match_list all_14242_genomes_contig_length_up.txt) Pub_cluster_partial_match_list

```
Pub_2946064051___Ga0505392_01
```
grep 'Pub_2946064051___Ga0505392_01' Pub_cluster_partial_match_list 
```
Pub_2946064051___Ga0505392_01
```
grep 'Pub_2946064051' all_14242_genomes_contig_length_up.txt
```
Pub_2946064051	Ga0505392_01	8957905
```
sed 's,Ga0505392_01,Pub_2946064051___Ga0505392_01,g' all_14242_genomes_contig_length_up.txt > all_14242_genomes_contig_length_up2.txt
grep -F -f Pub_cluster_partial_match_list all_14242_genomes_contig_length_up2.txt > all_14242_genomes_contig_length_up_partial_match_subset.txt # 31979
wc -l  Pub_cluster_partial_match_list #31979

# --------- now the contigID could match between cluster BGC file and contig length ---------
### the exact whole word match contig length file could be used directly to inner join with that of the BGC file (i.e., all_14242_genomes_contig_length_up_exact_word_match_subset.txt). However, for the partial matched cluster contigs, I need to generate a map file between Pub_cluster_partial_info.txt and the corresponding contigID within the contig length file 

for file in $(cat Pub_cluster_partial_match_list); do echo $file | tr '\n' '\t' >> all_14242_genomes_contig_length_up_partial_match_subset_up.txt; grep $file all_14242_genomes_contig_length_up_partial_match_subset.txt >> all_14242_genomes_contig_length_up_partial_match_subset_up.txt; done
wc -l all_14242_genomes_contig_length_up_partial_match_subset_up.txt # 31980
wc -l Pub_cluster_partial_match_list # 31979


###### ------- trouble shooting this one row difference --------

grep -v -w -F -f Pub_cluster_partial_match_list <(cut -f 1 all_14242_genomes_contig_length_up_partial_match_subset_up.txt | sort | uniq)
```
Pub_2511231008
```
grep -P '^Pub_2511231008\t' all_14242_genomes_contig_length_up_partial_match_subset_up.txt

```
Pub_2511231008	Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131Pseudomonassp.GM21:PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131	111952
```
grep 'Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131'  Pub_cluster_partial_match_list
grep 'Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131' all_14242_genomes_contig_length_up_partial_match_subset_up.txt

```
Pub_2511231008	Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131Pseudomonassp.GM21:PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131	111952
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131	Pub_2511231008	Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131Pseudomonassp.GM21:PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131	111952
```

grep -P -A 3 -B 3 'Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_1.1' Pub_cluster_partial_match_list # since the problemetic lines lack of its first column, and its former line is Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_1.1, I suspected that this may caused by some empty line within the Pub_cluster_partial_match_list, so I checked this. It turned out this is caused by the fact that when grep Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_1.1 it also catches Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131. due to the . wild character within the query

```
Pub_2511231007___PMI21_GM18_TGACCA_L004_R1_006_paired_trimmed_paired_contig_5.5
Pub_2511231007___PMI21_GM18_TGACCA_L004_R1_006_paired_trimmed_paired_contig_62.62
Pub_2511231007___PMI21_GM18_TGACCA_L004_R1_006_paired_trimmed_paired_contig_8.8
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_1.1
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_122.122
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_125.125
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_166.166
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_19.19
Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_20.20
```

#### ---- it turned out the Pub_2511231008___PMI22_GM21_ACAGTG_L004_R1_006_paired_trimmed_paired_contig_131.131 contigID has two line matches within all_14242_genomes_contig_length_up_partial_match_subset_up.txt. So I manually deleted the top one 

wc -l all_14242_genomes_contig_length_up_partial_match_subset_up.txt # 31979
wc -l Pub_cluster_partial_match_list # 31979

# ----------------  generate the contigID list that has length >=5Kb and then use this query to fish the final cluster infor file --------

cat <(awk '$4*1>=5000' all_14242_genomes_contig_length_up_partial_match_subset_up.txt| cut -f 1) <(awk '$3*1>=5000' all_14242_genomes_contig_length_up_exact_word_match_subset.txt | cut -f 2) > Iso_MAG_Pub_5kbplus_contigID_list # 28563 + 52377 = 80940
cat Iso_cluser_info_up2.txt MAG_cluster_info_up3.txt Pub_cluster_info_up3.txt > Iso_MAG_Pub_BGC_cluster_info_updated.txt
# 134489
grep -w -F -f Iso_MAG_Pub_5kbplus_contigID_list Iso_MAG_Pub_BGC_cluster_info_updated.txt  > Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt 
cut -f 3  Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt | sort | uniq | wc -l # 80940

####### ----- finally this BGC cluster file is subset to filter off those with contig length smaller than 5K #########

sed -i 's,Iso_Wt_WTC89_2___NODE,Iso_Wt_WTC89-2___NODE,g' Iso_MAG_Pub_BGC_cluster_5Kplus_info_updated.txt 
