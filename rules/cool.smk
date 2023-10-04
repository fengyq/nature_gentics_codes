
# generate min cool files.
rule bin_cooler:
    input:
        pair1=rules.pair.output,
    output:
        cool="results/cooler/{sample}_uniq.cool",
    log:
        "results/cooler/{sample}_uniq.cool.log"
    params:
        binpath= config["binpath"],
        chrom_sizes=config["chrom_sizes"],
        genome=config["assembly"]
    shell:
        """
        {params.binpath}bgzip -cd -@ 3 {input.pair1} |  {params.binpath}cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {params.genome} {params.chrom_sizes}:1000 - {output.cool} 2>{log}
        """


# generate multiresolution cool files.
rule zoom_cool:
    input:
        cool=rules.bin_cooler.output.cool,
    output:
        mcool="results/cooler/{sample}_uniq.mcool",
    log:
        "results/cooler/{sample}_uniq.mcool.log"
    params:
        binpath= config["binpath"],
        chrom_sizes=config["chrom_sizes"],
        res=config["resolutions"],
        balance=config["balance"]
    threads: 8
    shell:
        """
        {params.binpath}cooler zoomify --out {output.mcool} --resolutions {params.res} {params.balance} -p {threads} {input.cool} 2>{log}
        """


# merge cool files
# cooler merge A1.cool A2.cool A12.cool
# cooler zoomify --out A12.mcool --resolutions 1000,2000,5000,10000,50000,100000 --balance -p 16 A12.cool 2>{log}
# pairtools stats --merge ${stats} -o ${library_group}.$.stats
# merge multiple input stats files instead of calculating statistics of a .pairs/.pairsam file.
