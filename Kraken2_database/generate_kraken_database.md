## Build the Kraken database

### 1, Kraken microbiome database used for evaluate metagenomic bacteria ratio

##### a. build kraken2 default microbiome database
```
## library overview 
  1. GTDB bacteria and archaea QS>50 and low contam
  2. RefSeq viral,fungi, protozoa
  3. IMG_V
  4. FungiDB
  5. Ensemble_protist


## download Kraken2 default taxonomy
  DBNAME=$db/Kraken2_ncbi_taxdump_Nov23/Kraken_default
  kraken2-build --download-taxonomy --db $DBNAME

## download Kraken default microbiome libraries

  kraken2-build --download-library bacteria --db $DBNAME
  kraken2-build --download-library archaea --db $DBNAME
  kraken2-build --download-library protozoa --db $DBNAME
  kraken2-build --download-library fungi --db $DBNAME
  kraken2-build --download-library virus --db $DBNAME

## Build kraken2 default microbiome database

  mkdir log
  /mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --build --db $db/Kraken2_ncbi_taxdump_Nov23/Kraken_default --threads 96 >> log/Kraken2_default_build.log 2>&1
```
#### b. format GTDB genome fasta to match with that of Kraken2 requirement

```
## download GTDB genome and metadata

  wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz ./
  cat <(du -a GCA/ | grep 'fna.gz' | awk '{print $2}' |sed 's,/GCA,/\tGCA,g') <(du -a GCF/ | grep 'fna.gz' | awk '{print $2}' | sed 's,/GCF,/\tGCF,g') > genome_paths_FL_generated.tsv                                                    
  
  IFS=$'\n'
  for line in $(cat genome_paths_FL_generated.tsv)
  do
      output=$(echo $line | awk '{print $2}' | sed 's,.gz,,g')
      col2=$(echo $line | awk '{print $2}')
      col1=$(echo $line | awk '{print $1}')
      path=$(echo $col1$col2)
      echo $path
      gunzip -c $path > all_unzip_genome/$output
  done
  
  cat bac120_metadata_r207.tsv <(grep -v '^accession' ar53_metadata_r207.tsv) > bac120_and_ar53_metadata_r207.tsv
  cut -f -1 bac120_and_ar53_metadata_r207.tsv | sed 's,GB_,,g' | sed 's,RS_,,g' | sort | uniq > bac120_and_ar53_metadata_r207_accession_list

##  format taxonomy

	cat bac120_metadata_r207.tsv <(grep -v '^accession' ar53_metadata_r207.tsv) > bac120_and_ar53_metadata_r207.tsv
	cut -f 15 bac120_and_ar53_metadata_r207.tsv | sort | uniq | wc -l # 65704-1 rep species genomes 
	cut -f -1 bac120_and_ar53_metadata_r207.tsv | sed 's,GB_,,g' | sed 's,RS_,,g' | sort | uniq > bac120_and_ar53_metadata_r207_accession_list
	#### checked header of bacteria and archaea are the same
	
	cd all_unzip_genome 
	grep -w -F -f <(cat <(ls GCA_*.fna) <(ls GCF_*.fna) | sed 's,_genomic.fna,,g')  ../bac120_and_ar53_metadata_r207_accession_list | wc -l #65703
	cat <(ls GCA_*.fna) <(ls GCF_*.fna) | wc -l #65703, this is the total species genomes number of GTDBv207
	
	#### check if all ncbi taxid included in the meta file are included in the taxdum file in Kraken taxonomy file
	grep -v -w -F -f <(cut -f 1 ~/db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy/names.dmp  | sed 's, .*$,,g' | uniq) <(cut -f 1,78 bac120_and_ar53_metadata_r207.tsv | sed 's,^GB_,,g' | sed 's,^RS_,,g' | grep -v '^accession' | cut -f 2) # it turned out, lots of the taxid are not included in the ~/db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy/names.dmp.
	# subset genomes that has taxid matched with that in names.dmp,
	grep -w -F -f <(grep -w -F -f <(cut -f 1 ~/db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy/names.dmp  | sed 's, .*$,,g' | uniq) <(cut -f 1,78 bac120_and_ar53_metadata_r207.tsv | sed 's,^GB_,,g' | sed 's,^RS_,,g' | grep -v '^accession' | cut -f 2)) <(cut -f 1,78 bac120_and_ar53_metadata_r207.tsv | sed 's,^GB_,,g' | sed 's,^RS_,,g' | grep -v '^accession') > bac120_and_ar53_metadata_r207_matched_names_dmp_genomeID_vs_taxid_map.txt

	grep -w -F -f <(grep -v -w -F -f <(cut -f 1 ~/db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy/names.dmp  | sed 's, .*$,,g' | uniq) <(cut -f 1,78 bac120_and_ar53_metadata_r207.tsv | sed 's,^GB_,,g' | sed 's,^RS_,,g' | grep -v '^accession' | cut -f 2)) <(cut -f 1,78,79 bac120_and_ar53_metadata_r207.tsv | sed 's,^GB_,,g' | sed 's,^RS_,,g' | grep -v '^accession') > bac120_and_ar53_metadata_r207_not_match_names_dmp_meta_subset.txt
	
	paste -d  "\t"  <(cut -f 1 bac120_and_ar53_metadata_r207_not_match_names_dmp_meta_subset.txt)  <(cut -f 3 bac120_and_ar53_metadata_r207_not_match_names_dmp_meta_subset.txt | sed 's,.__,,g' | sed 's,;*$,,g'| sed 's,^.*;\(.*\)$,\1,g') > bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map.txt
	
	sed 's,\[Propionibacterium\] humerusii,Cutibacterium modestum,g' bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map.txt > bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map_up.txt
	
	paste -d "\t" <(cut -f 1 bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map_up.txt) <(taxonkit name2taxid --data-dir $db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy  <(cut -f 2 bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map_up.txt)) > temp.txt

	# after check for the output, there are several taxa got two taxid:
	
	#grep -w -F -f <(grep -w -F -f <(cut -f 2,3 temp.txt | sort | uniq | cut -f 1 | uniq -d) <(cut -f 2,3 temp.txt) | cut -f 1 | sort | uniq | cut -d ' ' -f 1 | sort | uniq) temp.txt  | cut -f 2,3 | sort | uniq | grep -v ' '
	#Bacillus	1386 # this is correct after grep from names.dmp file 
	#Bacillus	55087 - Bacillus <walking sticks> 
	#Bacteroidetes	200643 # this is class level 
	#Bacteroidetes	976 # this is correct, this is phylum level since the taxonomy in the GTDB meta file means phylum level
	#Gordonia	2053 # this is correct
	#Gordonia	79255 #Gordonia <eudicots>
	#Rhodococcus	1661425 #Rhodococcus <scale insects> 
	#Rhodococcus	1827 # this is correct

	paste -d "\t" <(cut -f 1 bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map_up.txt)  <(taxonkit name2taxid --data-dir $db/Kraken2_Feb11_2023/Mar23_dRep_7730_plusmore/taxonomy  <(cut -f 2 bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_map_up.txt) | grep -P -v 'Bacillus\t55087' | grep -P -v 'Bacteroidetes\t200643' | grep -P -v 'Gordonia\t79255' | grep -P -v 'Rhodococcus\t1661425') > bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_taxid_map_up.txt
	
	awk -F '\t' '{print $1 "\t" $3}' bac120_and_ar53_metadata_r207_not_match_names_dmp_genomeID_vs_taxnames_taxid_map_up.txt > bac120_and_ar53_metadata_r207_not_matched_names_dmp_genomeID_vs_taxid_map.txt
	cat bac120_and_ar53_metadata_r207_matched_names_dmp_genomeID_vs_taxid_map.txt  bac120_and_ar53_metadata_r207_not_matched_names_dmp_genomeID_vs_taxid_map.txt > bac120_and_ar53_metadata_r207_combined_genomeID_vs_taxid_map.txt 

	# ------- now the map file for each species genomes within all_unzip_genome folder has its taxid mapped ---
	cd all_unzip_genome 
	grep -w -F -f <(grep -w -F -f <(cat <(ls GCA_*.fna) <(ls GCF_*.fna) | sed 's,_genomic.fna,,g')  ../bac120_and_ar53_metadata_r207_accession_list ) ../bac120_and_ar53_metadata_r207_combined_genomeID_vs_taxid_map.txt > ../bac120_and_ar53_metadata_r207_combined_genomeID_vs_taxid_map_65703_species_sub.txt # subset map file to 65703 representative species set.
	
	cd $db/GTDB_database/gtdb_genomes_reps_r207
	
	IFS=$'\n'

	for file in $(cat bac120_and_ar53_metadata_r207_combined_genomeID_vs_taxid_map_65703_species_sub.txt)
	do
		echo $file
		genomeID=$(echo $file | cut -f 1)
		taxid=$(echo $file | cut -f 2)
		seqkit replace -p .+ -r "$genomeID"___Contig{nr} --nr-width 7 all_unzip_genome/"$genomeID"_genomic.fna > Kraken_formated_genomes/"$genomeID"_formated.fasta
		sed -i "s,${genomeID},${genomeID}|kraken:taxid|${taxid},g;s,___Contig, Contig,g" Kraken_formated_genomes/"$genomeID"_formated.fasta
	done
	
	## ----- cat all formated fasta file into one -----------
	cd $db/GTDB_database/gtdb_genomes_reps_r207/Kraken_formated_genomes # it is not possible to cat all files together list too long for ls and cat command
	cat GCF_*.fasta > All_GCF_species_rep_genome.fasta
	cat GCA_*.fasta > All_GCA_species_rep_genome.fasta
	cat All_GCF_species_rep_genome.fasta All_GCA_species_rep_genome.fasta > ../All_kraken_format_65703_species_rep_genomes.fasta
```

#### c. Build GTDB bacteria and archaea database with low contamination

```
## compare GTDB lowcontamination (QS>=50 and contamination<1%) with Kraken2 default Refseq bacteria and archaea genome

	cp -r $db/Kraken2_ncbi_taxdump_Nov23/library/* library/ # since GTDB database contains genomes from NCBI RefSeq genome, so, we curated GTDB genome in order to make it compatible with Kraken2 default database (RefSeq genomes)

#### RefSeq GenomeIDs included for Kraken build
	cat /mnt/m4/liufang/db/old_Kraken_before_Nov13/Kraken2_Feb11_2023/May25_krakenDefault_IRBCwithIMG_IMGVR3_FungiDB_EnsemblProtist/library/archaea/manifest.txt | sed 's,/,\t,g' |  cut -f 7 | sed 's,\(^..._...........\).*,\1,g' > Kraken_default_RefSeq_archaea_genomeID
	cat /mnt/m4/liufang/db/old_Kraken_before_Nov13/Kraken2_Feb11_2023/May25_krakenDefault_IRBCwithIMG_IMGVR3_FungiDB_EnsemblProtist/library/bacteria/manifest.txt | sed 's,/,\t,g' |  cut -f 7 | sed 's,\(^..._...........\).*,\1,g' > Kraken_default_RefSeq_bacteria_genomeID
	cat Kraken_default_RefSeq_archaea_genomeID Kraken_default_RefSeq_bacteria_genomeID > Kraken_default_RefSeq_bacteria_and_archaea_genomeID # 34158
	
	ln -s $db/Kraken2_ncbi_taxdump_Nov23/library/bacteria/assembly_summary.txt RefSeq_bacteria_assembly_summary_Kraken2_sep_27_2022.txt 
	grep -c -P 'Complete Genome|Chromosome' Kraken_default_bacteria_RefSeq_bacteria_assembly_summary.txt # 33663
	ln -s $db/Kraken2_ncbi_taxdump_Nov23/library/archaea/assembly_summary.txt RefSeq_archaea_assembly_summary_Kraken2_sep_27_2022.txt 
	grep -c -P 'Complete Genome|Chromosome' Kraken_default_archaea_RefSeq_archaea_assembly_summary.txt #495
	cat RefSeq_bacteria_assembly_summary_Kraken2_sep_27_2022.txt  <(grep -v '^#'  RefSeq_archaea_assembly_summary_Kraken2_sep_27_2022.txt) > RefSeq_bacteria_and_archaea_assembly_summary_Kraken2_sep_27_2022.txt
	cut -f 1 RefSeq_bacteria_and_archaea_assembly_summary_Kraken2_sep_27_2022.txt | grep -v '^#' > RefSeq_sep27_2022_all_bac_arc_assembly_genomeID
	grep -P 'Complete Genome|Chromosome' RefSeq_bacteria_and_archaea_assembly_summary_Kraken2_sep_27_2022.txt | cut -f 1 > RefSeq_sep27_2022_as_kraken2_reference_bac_arc_assembly_genomeID

## download the most updated version of bacteria and archaea assembly_summary.txt
	wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
	mv assembly_summary.txt  RefSeq_archaea_assembly_summary_Nov27_2023.txt
	grep -P -c 'Complete Genome|Choromosome' RefSeq_archaea_assembly_summary_Nov27_2023.txt #562
	wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
	mv assembly_summary.txt RefSeq_bacteria_assembly_summary_Nov27_2023.txt
	grep -P -c 'Complete Genome|Choromosome' RefSeq_bacteria_assembly_summary_Nov27_2023.txt #36322

## compared the overlapped genomes between Kraken genome and GTDB QS>50&LowContam

	#### GTDB repSpecies genomes 
	ln -s $db/Kraken2_gtdb_taxdump_Nov22/prep_genomes/GTDB/GTDB_35135_lowContam_NRgenomeID_list ./
	awk '$1==$16' GTDB_317542_all_genomes_bac120_and_ar53_metadata_r207.tsv | wc -l #65703
	awk '$17=="t"' GTDB_317542_all_genomes_bac120_and_ar53_metadata_r207.tsv  | wc -l #65703
	qwk '$17=="t"' GTDB_317542_all_genomes_bac120_and_ar53_metadata_r207.tsv  | grep -c '^GB_GCA_' # 42580
	awk '$17=="t"' GTDB_317542_all_genomes_bac120_and_ar53_metadata_r207.tsv  | grep -c '^RS_GCF_' # 23123

	#### GTDB genomes that used for generate Kraken database are the 35135 genomes that listed in GTDB_35135_lowContam_NRgenomeID_list
	grep -c '^GCA_' GTDB_35135_lowContam_NRgenomeID_list # 19500
	grep -c '^GCF_' GTDB_35135_lowContam_NRgenomeID_list #15635
	
	#### GTDB QS>50&LowContam Representative genomes are overlapped with Kraken bacteria and archaea genomes for a tiny portion.
	
	#grep -w -F -f GTDB_35135_lowContam_NRgenomeID_list RefSeq_bacteria_and_archaea_assembly_summary_Kraken2_sep_27_2022.txt | wc -l # 15189 of the genomes in GTDB are contained with Kraken2 default database
	wc -l Kraken_default_RefSeq_bacteria_and_archaea_genomeID #34158
	grep -c -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_bacteria_and_archaea_genomeID #3802
	ln -s $db/NCBI_Refseq/Mar5_2021/temp/wget.sh all_refseq_assembly_ftp_site_list
	grep -v -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_bacteria_and_archaea_genomeID > Kraken_default_RefSeq_bacteria_and_archaea_genomeID_not_coverred_by_GTDB_lowContam_RepGenome_list
	grep -v -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_bacteria_genomeID > Kraken_default_RefSeq_bacteria_genomeID_not_coverred_by_GTDB_lowContam_RepGenome_list #30077
	grep -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_bacteria_genomeID > Kraken_default_RefSeq_bacteria_genomeID_coverred_by_GTDB_lowContam_RepGenome_list
	grep -v -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_archaea_genomeID > Kraken_default_RefSeq_archaea_genomeID_not_coverred_by_GTDB_lowContam_RepGenome_list #279
	grep -w -F -f GTDB_35135_lowContam_NRgenomeID_list Kraken_default_RefSeq_archaea_genomeID > Kraken_default_RefSeq_archaea_genomeID_coverred_by_GTDB_lowContam_RepGenome_list 
	grep -v -w -F -f Kraken_default_RefSeq_bacteria_and_archaea_genomeID GTDB_35135_lowContam_NRgenomeID_list > GTDB_lowContam_not_coverred_by_Kraken2_default_GenomeID_list
	grep -v -w -F -f Kraken_default_RefSeq_bacteria_and_archaea_genomeID GTDB_35135_lowContam_NRgenomeID_list > GTDB_lowContam_not_coverred_by_Kraken2_default_GenomeID_list
```
## Format fasta and build library

```
## subset GTDB LowContam genome

	ln -s $db/Kraken2_gtdb_taxdump_Nov22/prep_genomes/GTDB/GTDB_35135_lowContam_NRgenomeID_list ./
	for file in $(cat GTDB_35135_lowContam_NRgenomeID_list); do echo $file; ln -s ../Kraken_formated_genomes/"$file"_formated.fasta ./; done
	cat *.fasta > GTDB_35135_lowContam_ncbi_formated_genomes_for_kraken.fasta
	
	ln -s $db/Kraken2_microbiome_Nov25/explore/compare_GTDB_lowContam_vs_RefSeq_Bac_Archaea/GTDB_lowContam_not_coverred_by_Kraken2_default_GenomeID_list ./
	for file in $(cat GTDB_lowContam_not_coverred_by_Kraken2_default_GenomeID_list); do echo $file; cp -s $db/GTDB_database/gtdb_genomes_reps_r207/Kraken_formated_QS50plus_lowContam/"$file"_formated.fasta ./; done
	cat *_formated.fasta > GTDB_lowContam_not_coverred_by_Kraken2_default_formated.fasta

## add library to Kraken default micorbiome database

	conda activate Kraken2.11
	cd $db/Kraken2_microbiome_Nov25
	/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --add-to-library $db/Kraken2_microbiome_Nov25/prep_genome/GTDB_ncbi_taxdump/GTDB_lowContam_not_coverred_by_Kraken2_default_formated.fasta --db $db/Kraken2_microbiome_Nov25 --threads 96 >> log/Kraken2_add_GTDB_lowContam_not_coverred_by_Krakendefault_to_library.log 2>&1
	cd $db/Kraken2_microbiome_Nov25/library
	mv added GTDB_not_cover_by_kraken_default
```
##  add IRVC to Kraken default microbiome database

```
## prepare genome and format curation
cd $db/CRVC_virus
ln -s $wd/IRBC_viromes/viral_sequence/result_231026/NR_viral_genomes/NR_CRVC_genomes.fna ./
sed 's,$, ,g' NR_CRVC_genomes.fna | sed 's, ,|kraken:taxid|10239 ,g' > NR_CRVC_genomes_formated.fna 

## add IRVC to the library
cd $db/Kraken2_microbiome_Nov25
/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --add-to-library $db/CRVC_virus/NR_CRVC_genomes_formated.fna --db $db/Kraken2_microbiome_Nov25 --threads 96 >> log/Kraken2_add_CRVC.log 2>&1
cd $db/Kraken2_microbiome_Nov25/library
mv added CRVC
```

## add FungiDB genome reference into Kraken microbiome

```
## prepare genome


## build the library and add to Kraken microbiome
cd $db/Kraken2_microbiome_Nov25
/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --add-to-library $db/FungiDB/genomes/combine/All_Fungi_genome.fasta --db  $db/Kraken2_Feb11_2023 --threads 96 >> log/Kraken2_add_fungi_library.log 2>&1                                                                                                           
mv added/ FungiDB
```

## add ENSEMBLE protist gneome reference to Kraken microbiome

```
## prepare genome and fomat the fasta file

## build the library
cd $db/Kraken2_Feb11_2023
/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --add-to-library $db/EnSEMBL_protist/formated_genomes/uniq_header_fasta/All_EnSEMBL_protist_genome.fasta --db  $db/Kraken2_Feb11_2023 --threads 64 >> log/Kraken2_add_protist_library.log 2>&1

```


## build Kraken2 defaul library

```
conda activate Kraken2.11
/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --build --fast-build --db $db/Kraken2_microbiome_Nov25 --threads 132 >> log/Kraken2_microbiome_fast_build.log 2>&1

/mnt/m1/liufang/anaconda3/envs/Kraken2.11/bin/kraken2-build --build --db $db/Kraken2_microbiome_Nov25 --threads 132 >> log/Kraken2_microbiome_build.log 2>&1
```
### 1. Kraken2 default database - classification improvement by CRBD
```
#TOC
#1. Kraken2 default Sep27_2022
$db/Kraken2_ncbi_taxdump_Nov23/Kraken_default
#2. Kraken2 default + PubCrop
$db/Kraken2_ncbi_taxdump_Nov23/Kraken_default_plus_PubCrop
#3. Kraken2 defalut + PubCrop + CRBC
$db/Kraken2_ncbi_taxdump_Nov23/Kraken_default_plus_PubCrop_CRBC
#4. GTDB_PubCrop R207 lowContam MQ
$db/Kraken2_gtdb_taxdump_Nov22/GTDB_default
#5. GTDB_PubCrop R207 + PubCrop
$db/Kraken2_gtdb_taxdump_Nov22/GTDB_plus_PubCrop
#6. GTDB_PubCrop R207 + PubCrop + CRBC
$db/Kraken2_gtdb_taxdump_Nov22/GTDB_plus_PubCrop_CRBC_RefSeq



$db/Kraken2_phage_RepSpecies
$db/Kraken2_microbiome_exclude_virus
$db/Kraken2_microbiome_Nov25
```

### 2. Kraken2 + CRBD + Fungi + virus

### 3. Kraken2 + CRBD + Fungi + virus + GTDB

### 4. Kraken2 + GTDB only
