##### generate_rulegraph #####
rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        "results/rulegraph.png"
    params:
        binpath= config["binpath"],        
    shell:
        """
        {params.binpath}snakemake --rulegraph --configfile config.yml | dot -Tpng > {output}
        """
