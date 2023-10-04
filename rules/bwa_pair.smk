# fastp : trim and filter out bad reads
# bwa_mem : mapping reads using bwa;

rule bwa_map:
    input:
        fastq1=rules.split_linker.output.n5,
        fastq2=rules.split_linker.output.n6
    output:
        json="results/bwa_map/{sample}_split.fastp.json",
        html="results/bwa_map/{sample}_split.fastp.html",
        bam="results/bwa_map/{sample}.F2316.bam"
    params:
        fastp_options=config["fastp_options"],
        bwa_options=config["bwa_options"],
        sam_options=config["samtools_options"],
        binpath= config["binpath"],
        bwa_index=config["bwa_index"]
    shell:
        """
        {params.binpath}fastp {params.fastp_options} --in1 {input.fastq1} --in2 {input.fastq2} --json {output.json} --html {output.html} --stdout  | {params.binpath}bwa mem {params.bwa_options} {params.bwa_index} - | {params.binpath}samtools view {params.sam_options}  > {output.bam}
        """

# statistic of bam
rule samstats:
    input:
        rules.bwa_map.output.bam
    output:
        "results/bwa_map/{sample}.F2316.bam.stat"
    params:
        binpath= config["binpath"],
    shell:
        """
        {params.binpath}samtools stat  {input} > {output}
        """


# make bedpe files for interleaved paired-bam file
rule bedpe:
    input:
        rules.bwa_map.output.bam
    output:
        "results/bwa_map/{sample}.F2316.bedpe.gz"
    params:
        binpath= config["binpath"],
    shell:
        """
        {params.binpath}bedtools bamtobed -bedpe -i  {input} |awk '{{if ($0 !~ /chrEBV|random|chrUn_/) print $0}}'  | bgzip > {output}
        """

# make sorted pair files in the lexicographic order
rule pair:
    input:
        rules.bedpe.output
    output:
        pair="results/pair/{sample}.uniq.pairs.gz",
        stat="results/pair/{sample}.pairs.stats"
    params:
        binpath= config["binpath"],
        header= config["pair_header"],
        tmpdir= "`mktemp -d`",
        sort_options= config["sort_options"]
    shell:
        """
        # use middle position of bedPE fragments as pairs
        zcat {input} | awk '{{print $7,$1,int(($2+$3)/2),$4,int(($5+$6)/2),$9,$10,"UU"}}' OFS="\\t" | cat {params.header} - | \
        {params.binpath}pairtools sort --tmpdir {params.tmpdir} {params.sort_options} | {params.binpath}pairtools dedup --max-mismatch 1 --output {output.pair}  --output-stats {output.stat} 
        """

# index pair file
rule pair_index:
    input:
        rules.pair.output.pair
    output:
        "results/pair/{sample}.uniq.pairs.gz.px2"
    log:
        "results/pair/{sample}.uniq.pairs.px2.log"
    params:
        binpath= config["binpath"]
    shell:
        """
        {params.binpath}pairix -p pairs {input} 2>{log}
        """
