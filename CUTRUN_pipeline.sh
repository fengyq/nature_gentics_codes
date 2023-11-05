# CUT and RUN pipeline
# @Author: feng_yuanqing

#----------------------------------------------
# 1.extract fastq from bcl
#----------------------------------------------
bcl2fastq --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  --no-lane-splitting \
  --barcode-mismatches 1 \
  -r 6 -w 6 \
  -R /project/AHNF3WBGXN \
  --output-dir=rawfastq2 --sample-sheet=/cutrun/samplesheet.csv


#----------------------------------------------
# 2. fastq qc
#----------------------------------------------
cd /cutrun/rawfastq0/
for x in $(ls *gz); do
  echo $x
  fastqc -q -t 4 --noextract -f fastq --outdir fastqc $x
done

#----------------------------------------------
# 3. filter fastq
#----------------------------------------------
cd /cutrun/rawfastq0/fastp
for x in $(cat files.tsv); do
  echo $x
  fastp --thread 6 --overrepresentation_analysis --detect_adapter_for_pe --in1 /cutrun/rawfastq0/${x}1_001.fastq.gz --in2 /cutrun/rawfastq0/${x}2_001.fastq.gz --out1 ${x}1_fastp.fq.gz --out2 ${x}2_fastp.fq.gz --length_required 20 --html ${x}.fastp.html 2>> ${x}.fastp.log 1>&2
done

#----------------------------------------------
# 4. cut and run qc: k-MetaSTAT nucleosome barcodes; H3K4Me3 Vs IgG
#----------------------------------------------
# template loop begin ##
echo "R1.fq 250bp 22nt barcodes = 8 counts"
for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA; do
  grep -c $barcode <(zcat mnt1-cut-H3K4Me3_R1_fastp.fq.gz)
done

echo "R2.fq 250bp 22nt barcodes = 8 counts"
for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA; do
  grep -c $barcode <(zcat mnt1-cut-IgG_S3_R1_fastp.fq.gz)
done
# template loop end ##

#----------------------------------------------
# 5. mapping using bowtie2
#----------------------------------------------
# --dovetail: concordant when mates extend past each other
# mapp reads to HG38
for x in $(cat files); do
  echo $x
  bowtie2 --threads 10 --dovetail --very-sensitive-local -I 10 -X 700 -x bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -1 /cutrun/rawfastq0/fastp/${x}_R1_fastp.fq.gz -2 /cutrun/rawfastq0/fastp/${x}_R2_fastp.fq.gz | samtools view -u - | samtools sort -@ 8 > ${x}.sorted.bam 2> ${x}.bowtie.log
  samtools index ${x}.sorted.bam
done
# mapp reads to ecoli; calculate scale factors
for x in $(cat files); do
  echo $x
  bowtie2 --threads 10 --dovetail --very-sensitive-local -I 10 -X 700 -x Ecoli_MG1655 -1 /cutrun/rawfastq0/fastp/${x}_R1_fastp.fq.gz -2 /cutrun/rawfastq0/fastp/${x}_R2_fastp.fq.gz | samtools view -u -@ 10 -h -b -q 20 -F 1804 -f 2 - | samtools sort -@ 8 > Ecoli_${x}.sorted.bam
  samtools index Ecoli_${x}.sorted.bam
done

#----------------------------------------------
# 6. Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
#----------------------------------------------
cd /cutrun/mnt1_Cutrun
for x in $(cat files); do
  echo $x
  java -Xmx10G -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar gatk4-4.1.7.0-0/gatk-package-4.1.7.0-local.jar MarkDuplicates --QUIET true --INPUT ${x}.sorted.bam --OUTPUT ${x}.dup_marked.bam --METRICS_FILE ${x}.dup_marked.metrics --REMOVE_DUPLICATES false --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
  samtools view -@ 10 -h -b -q 20 -F 1804 -f 2 ${x}.dup_marked.bam > ${x}.rmdup.bam
  samtools index -@ 10 ${x}.rmdup.bam
done


#----------------------------------------------
# 7. generate bigwig and bedGraph for visulization
#----------------------------------------------
bamCoverage -b mnt1-cut-xx.rmdup.bam -o mnt1-cut-xx.rmdup.bw -p 8 --scaleFactor 0.438 --normalizeUsing RPKM --binSize 10 --blackListFileName ~/genome/ATAC_blacklist/hg38-blacklist.v2.bed --extendReads


#----------------------------------------------
# 8. call peaks and calculate FRiP; overall SEAR is better for TF cut&run; use macs2 for broad peak without control.
#----------------------------------------------
# call broad peak
macs2 callpeak -t mnt1-cut-H3K27ac.rmdup.bam -g hs --broad --broad-cutoff 0.05 -f BAMPE -n mnt1-cut-H3K27ac --outdir mnt1-cut-H3K27ac --keep-dup all 2> mnt1-cut-H3K27ac.macs2.log

# call narrow peak
# https://github.com/FredHutch/SEACR;Since we have normalized fragment counts with the E. coli read count, we set the normalization option of SEACR to “non”. Otherwise, use “norm”
bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
bash SEACR_1.3.sh ../bedGraph/mnt1-cut-H3K4Me3.scaled.bedGraph ../bedGraph/mnt1-cut-IgG_S3.scaled.bedGraph non relaxed mnt1-cut-H3K4Me3_SEACR
bash SEACR_1.3.sh ../bedGraph/mnt1-cut-H3K4Me3.scaled.bedGraph ../bedGraph/mnt1-cut-IgG_S3.scaled.bedGraph non stringent mnt1-cut-H3K4Me3_SEACR
awk -v OFS="\t" '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1],region[1],region[2]}' mnt1-cut-H3K4Me3_SEACR.stringent.bed > mnt1-cut-H3K4Me3_SEACR.stringent.center

## Make a custom "SAF" file which featureCounts needs:
awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2, $3, "+"}' mnt1-cut-H3K4Me3_SEACR.relaxed.bed > mnt1-cut-H3K4Me3.saf
## run featureCounts (add -T for multithreading)
featureCounts -T 10 -p -a mnt1-cut-H3K4Me3.saf -F SAF -o mnt1-cut-H3K4Me3.saf.out ../rmdup_bam/mnt1-cut-H3K4Me3.rmdup.bam

#----------------------------------------------
# 9. TSS enrichment Analysis
#----------------------------------------------
for x in $(cat files); do
  echo $x
  fast_deep_plot2.sh -a /cutrun/mnt1_Cutrun/bigwig/$x.rmdup.bw
  mv ./deeptools-plot/plotHeatmap.matrix.gz ./deeptools-plot/$x.matrix.gz
done