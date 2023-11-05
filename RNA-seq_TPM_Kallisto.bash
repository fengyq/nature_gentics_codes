#Bulk RNA_seq analysis workflow
## step 0. check reads queality by fastqc, trim adapter if necessary.
## step 1. Kallisto. Mapping in MSI login system.
#1) index. get mRNA and index.
#2) mapping. make directory and run kallisto script, get abundance result file. 
##step 2. Tximport. Ensemble txid to gene name in Rstudio.
#1) index. get gene to refMrna (in MSI login system).
#2) tx2gene. tximport from kallisto, ensembel txid to gene name in Rstudio.
##step 3. DEBrowser. DEG analysis in Rstudio and website.

##------------------
1. Fastqc  :run_fastqc.sh 
##------------------
# fastqc 
cd RNA_Seq/CRISPRi
for x in $(ls *_001.fastq.gz); do echo "starting fastqc : ${x}"; fastqc --noextract --format fastq -t 10 -o RNA_Seq/fastqc $x ; done 

# multiqc report
conda activate hic
cd RNA_Seq/fastqc
multiqc ./fastqc --filename fastqc_output.report

##------------------
2. trim adapter and filter reads
##------------------
# fastp_filter.sh
cd RNA_Seq/CRISPRi
for x in $(ls *_001.fastq.gz | sed 's/_R[1-2]_001.fastq.gz//g' |uniq ); do echo "sample: ${x}"; fastp --thread 18 --in1 ${x}_R1_001.fastq.gz --in2 ${x}_R2_001.fastq.gz --out1 data/RNA_Seq/fastp_filtered/${x}_R1.trimed.fastq.gz --out2 data/RNA_Seq/fastp_filtered/${x}_R2.trimed.fastq.gz  --json data/RNA_Seq/fastp_filtered/${x}.json  --html data/RNA_Seq/fastp_filtered/${x}.html  2>>data/RNA_Seq/fastp_filtered/${x}_fastp.log 1>&2 ; done 
# bsub -n 12 -M 20000 -R "span [hosts=1]" -e myjob.err -o myjob.out -J fastp_trim sh fastp_filter.sh 


##------------------
3. kallisto pipiline
##------------------
## get mRNA and index
wget hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
kallisto index -i refMrna.fa.idx refMrna.fa.gz
##gene to refMrna
wget hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat <refGene.txt.gz |cut -f 2,13 | sort| uniq > tx2gene

##------ kallisto_run.sh------
cd RNA_Seq/fastp_filtered
for x in $(ls *fastq.gz | sed 's/_R[1-2].trimed.fastq.gz//g' |uniq ); 
do echo "sample: ${x}"; mkdir -p RNA_Seq/kallisto_tpm/${x};
kallisto quant -i genome/hg38_mRNA/refMrna.fa.idx -o RNA_Seq/kallisto_tpm/${x} -t 16 ./${x}_R1.trimed.fastq.gz ./${x}_R2.trimed.fastq.gz;
done 
# bsub -n 15 -M 20000 -R "span [hosts=1]" -e myjob.err -o myjob.out -J kallisto_quant sh kallisto_run.sh 



##-----------------------
Tximport pipiline
##-----------------------
##  txid to gene name 
Rscript RNA_Seq/kallisto_tpm/kallisto_tximport.R 
## or  Run the lines below in R or RStudio.
## tximport from kallisto
## convert txid to gene name
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Import_transcript-level_estimates
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("rhdf5")

# read in gene-transript table
library(tximport)
dir="~/genome/hg38_mRNA/"
tx2gene <- read.table(file.path(dir, "tx2gene"), header = TRUE)
head(tx2gene)

# read in sample name
dir2="RNA_Seq/kallisto_tpm"
samples <- read.table(file.path(dir2, "samples.txt"), header = TRUE)
files <- file.path(dir2, samples$samples, "abundance.tsv")
names(files) <- samples$samples

# summarize gene-level counts
library(rhdf5)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene) ## read in tsv files
head(txi$counts) ## counts at gene level
names(txi)

# format count table and write
library(tibble)
G_count<- round(txi$counts,2)
FINAL_count<-as_tibble(G_count,rownames = "id")

write.table(FINAL_count, file= paste0("kallisto_gene.rawcount"), sep = "\t", row.names = F, col.names = T,quote = F)

G_tpm<- round(txi$abundance,2)
FINAL_tpm<-as_tibble(G_tpm,rownames = "id")

write.table(FINAL_tpm, file= paste0("kallisto_gene.tpm"), sep = "\t", row.names = F, col.names = T,quote = F)

##-----------------------
DEBrowser pipiline
##-----------------------

##Install debrowser
# Installation instructions:
# 1. Install DEBrowser and its dependencies by running the lines below in R or RStudio.

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("debrowser")

# 2. Load the library

library(debrowser)

# 3. Start DEBrowser

startDEBrowser()

