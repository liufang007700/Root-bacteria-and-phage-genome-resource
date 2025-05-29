## a. run antismash6.1.0
##### This section is done by Ray
```bash
screen -r antismash2
source /mnt/m2/dairui/anaconda3/bin/activate
conda activate antismash

DIR=/mnt/m2/dairui/project/binning/isolate/
# --minlength default 1000
# --allow-long-headers 没有这个的话会报错
mkdir -p $DIR/08_BGCs/Wheat/sep_genome

time parallel -j 2 --xapply \
antismash \
  -c 12 \
  --allow-long-headers \
  --asf \
  --cb-general \
  --cb-knownclusters \
  --cb-subclusters \
  --pfam2go \
  --output-dir $DIR/08_BGCs/Wheat/sep_genome/{1} \
  --genefinding-gff3 $DIR/04_prodigal/Wheat/sep_genome/{1}.gff \
  $DIR/01_genome/Wheat/{1}.fa \
  ::: `cat $DIR/01_genome/Iso_Wt_genomeID`
 # 1599294.55s user 463881.16s system 519% cpu 110:12:52.80 total
```
