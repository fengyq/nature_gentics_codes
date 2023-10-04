# @Author: Yuanqing Feng <fengyq>
# @Date:   Monday, July 12th 2021, 8:09:09 pm
# @Email:  fengyq2018@gmail.com
# @Filename: fastqc.smk
# @Last modified by:   fengyq
# @Last modified time: Monday, July 12th 2021, 8:12:47 pm
# @Copyright: fengyq2018@gmail.com



# Run FastQC on a FASTQ file
rule fastqc:
    input:
        "raw_data/{id}.fq.gz"
    output:
        "results/fastqc/{id}_fastqc.html",
        "results/fastqc/{id}_fastqc.zip"
    params:
        binpath= config["binpath"],
    shell:
        """
        {params.binpath}fastqc {input} -q -o ./results/fastqc
        """

# Aggregate all FastQC reports into a MultiQC report.
rule multiqc:
    input:
        expand("results/fastqc/{id}_fastqc.zip", id = config["fastq_ids"]),
        expand("results/bwa_map/{sample}_split.fastp.json", sample=config["sample_ids"]),
        expand("results/bwa_map/{sample}.F2316.bam.stat", sample=config["sample_ids"]),
        # using id as wildcard to match all samples; iinput end with .fq.gz
    output:
        html = "results/multiqc/multiqc.html",
        log = "results/multiqc/multiqc.log"
    params:
        binpath= config["binpath"],
    shell:
        """
        # Run multiQC and keep the html report
        {params.binpath}multiqc  {input} -o ./results/multiqc/ 2> {output.log}
        """
