# @Author: Yuanqing Feng <fengyq>
# @Date:   Monday, July 12th 2021, 8:09:09 pm
# @Email:  fengyq2018@gmail.com
# @Filename: split_linker.smk
# @Last modified by:   fengyq
# @Last modified time: Thursday, September 30th 2021, 7:19:42 pm
# @Copyright: fengyq2018@gmail.com


# Run Flash to merge two paired fastqs
rule Flash_merge:
    input:
        r1="raw_data/{sample}_R1.fq.gz",
        r2="raw_data/{sample}_R2.fq.gz",
    output:
        f0="results/Flash/{sample}_flash.extendedFrags.fastq", # infer wildcard from rule_all
        f1="results/Flash/{sample}_flash.notCombined_1.fastq",
        f2="results/Flash/{sample}_flash.notCombined_2.fastq",
    params:
        flash = config["flash"],
        binpath= config["binpath"],
    shell:
        """
        {params.binpath}flash {input.r1} {input.r2} {params.flash} -o {wildcards.sample}_flash -d results/Flash 1>results/Flash/{wildcards.sample}_flash.log
        """

# Run Flash to merge two paired fastqs
rule split_linker:
    input:
        f0=rules.Flash_merge.output.f0,
        f1=rules.Flash_merge.output.f1,
        f2=rules.Flash_merge.output.f2,
    output:
        n0="results/Split/{sample}_nolinker.fastq",
        n1="results/Split/{sample}_new1.fastq", # splitted read1 from flash merged read
        n2="results/Split/{sample}_new2.fastq", # splitted read2 from flash merged read
        n3="results/Split/{sample}_new3.fastq", #  read1 from flash unmerged read
        n4="results/Split/{sample}_new4.fastq", #  read2 from flash unmerged read
        n5="results/Split/{sample}_split1.fq.gz", #  read1
        n6="results/Split/{sample}_split2.fq.gz", #  read2
        stat1="results/Split/{sample}_linker.stat.mergedSE.tsv",
        stat2="results/Split/{sample}_linker.stat.PE.tsv",
        log1="results/Split/{sample}_linker_split.mergedSE.log",
        log2="results/Split/{sample}_linker_split.PE.log",
    shell:
        """
        # split merged fastq
        python ./code/split_linker_mergedSE.py -M 18 -S {output.stat1} -L ".CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG." -r1 {input.f0} -R0 {output.n0} -R1 {output.n1} -R2 {output.n2} 2>&1 | tee -a  {output.log1}
        # split the two unmerged fastqs
        python ./code/split_linker_PEv2.py -M 18 -S {output.stat2} -L ".CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG." -r1 {input.f1} -r2 {input.f2} -R1  {output.n3} -R2 {output.n4}  2>&1 | tee -a {output.log2}
        #merge all 4 fastqs
        cat {output.n1} {output.n3} |bgzip > {output.n5}
        cat {output.n2} {output.n4} |bgzip > {output.n6}
        """
