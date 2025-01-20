## Dissect the PGPR potential among CRBD

### 1. Conduct genome function annotation using diamond

#### a. Build diamond database
```

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
