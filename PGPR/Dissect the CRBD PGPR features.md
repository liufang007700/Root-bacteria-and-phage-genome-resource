## Dissect the PGPR potential among CRBD

### 1. Conduct genome function annotation using diamond

#### a. Build diamond database
```
## first build a kegg protein database including all eukaryotes, prokaryotes and virus, please be noted the kegg database we used is the April 1st 2021 version
cd kegg/genes/fasta
gunzip eukaryotes.pep.gz
gunzip prokaryotes.pep.gz
gunzip T40000.pep.gz
cat eukaryotes.pep prokaryotes.pep T40000.pep > KEGG_eukaryotes_prokaryotes_virus_protein.faa

## extract peptide that encode PGPR functions, for detailed PGPR function defination please refer to the supplenmental table

### prepare all PGPR peptide sequences

grep -w -F -f PGPR_KO_list ~/db/ftp_download_KEGG_04_01_2021/ko_gene_map.txt  | cut -f 1 > All_PGPR_KeggGeneID_list #153498 KeggGeneID
seqkit grep -f All_PGPR_KeggGeneID_list  KEGG_eukaryotes_prokaryotes_virus_protein.faa > All_PGPR_all_organisms_KeggGene_peptide.faa # 

### add CPS and KS protein into above (note: cause lots of bacteria are able to synthesize GA, however, kegg database did not include as much as GA synthesis related proteins, so we manually download and add into the PGPR protein reference)
    > CPS accession ID:NP_768789, NP_106893, NP_443949, NP_659791, WP_003466962, AEQ94336, WP_020322919
    > KS accession ID: NP_768790,NP_106894,NP_443948,NP_659792,WP_003466963,AEQ94335
conda activate entrez_direct 
esearch -db protein -query "$(cat KS_accession_list | tr '\n' ' ')" | efetch -format fasta > KS_protein_sequences.fasta
esearch -db protein -query "$(cat CPS_accession_list | tr '\n' ' ')" | efetch -format fasta > CPS_protein_sequences.fasta

### Add CPS and KS as prefix in the header line for both faa file and then combine together to construct diamond db
cat CPS_protein_sequences.fasta  KS_protein_sequences.fasta > CPS_KS_protein_seq.faa

## make diamond database

mkdir log
time \
diamond makedb \
--in All_PGPR_all_organisms_KeggGene_peptide.faa \
--db /mnt/m1/liufang/db/PGPR_DIAMOND/PGPR_dmnd \
--threads 64 \
>> log/diamond_makedb.log 2>&1

```
#### b. Conduct blast

```
mkdir log
mkdir diamond_out
time \
    /mnt/m2/dairui/anaconda3/envs/diamond2.0.15/bin/diamond blastp \
    --query CRBD_all_9772_prodigal_combine.faa \
    --db /mnt/m1/liufang/db/PGPR_DIAMOND_2023/All_Kegg_PGPR_plus_pub_CPS_KS_dmnd.dmnd \
    --threads 96 \
    --out diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6 \
    --outfmt 6 \
    --log \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --sensitive \
    --block-size 4 \
    --index-chunks 1 \
    >> log/CRBD_all_9772_genomes_faa_vs_PGPR_dmndblastp.log 2>&1

## Since the manually added KO 
cat ko_gene_map.txt Bacteria_CPS_KS_KO_gene_map.txt > ko_gene_map_add_pub_CPS_KS.txt
paste -d "\t" <(sed 's,.*___Pub,Pub,g' diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6 | cut -f 1 | sed 's,___,\t,g' | cut -f 1) <(sed 's,.*___Pub,Pub,g' diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6) > CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot_up.fm6
grep -w -v -F -f Phosphatase_SCL_results_link/KO_list_used_for_SCL_prediction CRBD_all_9772_genomes_PGPR_Dec20_evalue5_id50plus_KO_count.txt > CRBD_all_9772_genomes_PGPR_Dec20_exclude_phosphorus_SCL_evalue5_id50plus_KO_count.txt
grep -w -F -f Phosphatase_SCL_results_link/KO_list_used_for_SCL_prediction CRBD_all_9772_genomes_PGPR_Dec20_evalue5_id50plus_KO_count.txt  | cut -f 1 # those below 9 KOs are removed since their SCL KO listed were updated based on psortb results
```
