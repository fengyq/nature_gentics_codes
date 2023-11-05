# @Author: feng_yuanqing
# ATAC-seq pipeline 

#----------------------------------------------
# 1. fastq qc
#----------------------------------------------
for x in $(ls *gz); do
  echo $x
  fastqc -q -t 4 --noextract -f fastq --outdir fastqc $x
done


#----------------------------------------------
# 2. filter fastq
#----------------------------------------------

cd atac/rawfastq0/fastp
for x in $(cat files.tsv); do
  echo $x
  fastp --thread 6 --overrepresentation_analysis --detect_adapter_for_pe --in1 atac/rawfastq0/${x}1_001.fastq.gz --in2 atac/rawfastq0/${x}2_001.fastq.gz --out1 ${x}1_fastp.fq.gz --out2 ${x}2_fastp.fq.gz --length_required 20 --html ${x}.fastp.html 2>> ${x}.fastp.log 1>&2
done


#----------------------------------------------
# 3. ATAC mapping using bowtie2
#----------------------------------------------
cd atac/bowtie_map
for x in $(cat atac.files); do
  echo $x
  bowtie2 --threads 10 --very-sensitive --maxins 1000 -x bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -1 atac/rawfastq0/fastp/${x}_R1_fastp.fq.gz -2 atac/rawfastq0/fastp/${x}_R2_fastp.fq.gz | samtools view -u - | samtools sort -@ 4 > ${x}.sorted.bam 2> ${x}.bowtie.log
  samtools index ${x}.sorted.bam
done

#----------------------------------------------
# 4. mark duplicate reads
#----------------------------------------------
cd atac/bowtie_map
conda activate gatk4
java -Xmx10G -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar gatk-package-4.1.7.0-local.jar MarkDuplicates --QUIET true --INPUT ATAC.sorted.bam --OUTPUT ATAC.dup_marked.bam --METRICS_FILE ATAC.dup_marked.metrics --REMOVE_DUPLICATES false --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT

#----------------------------------------------
# 5. Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
#----------------------------------------------
samtools view -@ 6 -h -b -q 20 -F 1804 -f 2 ATAC.dup_marked.bam > ATAC.rmdup.bam
samtools index -@ 6 ATAC.rmdup.bam

#----------------------------------------------
# 6. generate bigwig (RPKM normalizztion)for visulization
#----------------------------------------------
bamCoverage -b ATAC.rmdup.bam -o ATAC.rmdup.bw -p 8 --normalizeUsing RPKM --binSize 10 --blackListFileName ~/genome/ATAC_blacklist/hg38-blacklist.v2.bed --extendReads
bamCoverage -b ATAC.rmdup.bam -o ATAC.rmdup.bw -p 8 RPGC --effectiveGenomeSize 2913022398 --binSize 10 --blackListFileName ~/genome/ATAC_blacklist/hg38-blacklist.v2.bed --extendReads

#----------------------------------------------
# 7. call peaks and calculate FRiP (fraction of reads in peaks )
#----------------------------------------------
# macs2 narrow peak calling
macs2 callpeak -t ATAC.rmdup.bam -g hs -f BAMPE -n ATAC --outdir macs2 -q 0.01 --keep-dup all 2>&1 1> macs2.log

## Make a custom "SAF" file which featureCounts needs:
awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' ATAC_peaks.narrowPeak > ATAC_peaks.saf
## run featureCounts (add -T for multithreading)
featureCounts -T 10 -p -a ATAC_peaks.saf -F SAF -o ATAC_peaks.saf.txt ../ATAC.rmdup.bam

#----------------------------------------------
# 8. peak enrichment in TSS
#----------------------------------------------
gtftools -g gencode.v42.gene.bed gencode.v42.basic.annotation.gtf.ensembl     # Make gene bed file
fast_deep_plot2.sh -a atac/ATAC/mnt-ATAC.rmdup.bw  # enrichment plot

#plot correlation heatmap
multiBigwigSummary bins -b atac/ATAC/mnt-ATAC.rmdup.bw atac/wm88-ATAC/wm88-ATAC.rmdup.bw -o mnt-wm88-ATAC.corr.npz
plotCorrelation -in mnt-wm88-ATAC.corr.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of bigwig Scores" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mnt-wm88-ATAC_SpearmanCorr.pdf --outFileCorMatrix mnt-wm88-ATAC_SpearmanCorr.tab

