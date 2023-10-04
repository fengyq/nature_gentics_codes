# @Author: Yuanqing Feng <fengyq>
# @Date:   2021-07-06T10:08:27-04:00
# @Email:  fengyq2018@gmail.com
# @Last modified by:   fengyq
# @Last modified time: Thursday, September 30th 2021, 7:19:37 pm
# @Copyright: fengyq2018@gmail.com

# this is a template from https://github.com/NBISweden/workshop-reproducible-research/tree/main/docker

from snakemake.utils import min_version
min_version("5.3.0")

configfile: "config.yml"

# Collect the main outputs of the workflow.
rule all:
    input:
        "results/multiqc/multiqc.log",
        "results/rulegraph.png",
        expand("results/Flash/{sample}_flash.extendedFrags.fastq", sample=config["sample_ids"]), # define a wildcard name : sample
        expand("results/Split/{sample}_linker.stat.PE.tsv", sample=config["sample_ids"]), # define a wildcard name : sample
        expand("results/bwa_map/{sample}_split.fastp.json", sample=config["sample_ids"]),
        expand("results/bwa_map/{sample}.F2316.bam", sample=config["sample_ids"]),
        expand("results/bwa_map/{sample}.F2316.bam.stat", sample=config["sample_ids"]),
        expand("results/bwa_map/{sample}.F2316.bedpe.gz", sample=config["sample_ids"]),
        expand("results/pair/{sample}.uniq.pairs.gz",sample=config["sample_ids"]),
        expand("results/pair/{sample}.uniq.pairs.gz.px2",sample=config["sample_ids"]),
        expand("results/pair/{sample}.pairs.stats", sample=config["sample_ids"]),
        expand("results/cooler/{sample}_uniq.cool", sample=config["sample_ids"]),
        expand("results/cooler/{sample}_uniq.mcool", sample=config["sample_ids"]),

##### load rules #####
include: "rules/fastqc.smk",
include: "rules/general.smk",
include: "rules/split_linker.smk",
include: "rules/bwa_pair.smk",
include: "rules/cool.smk",


#### clean temp fastq results
rule clean:
    shell:
        """
        rm -f results/Split/*fastq; rm -f results/Flash/*fastq
        """
# you can run it by passing the rule name to snakemake:
# snakemake clean
