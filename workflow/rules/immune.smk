rule immunedeconv:
    input:
        "results/counts/all_tpm.tsv",
    output:
        mcp="results/immune/mcp.pdf",
        quantiseq="results/immune/quantiseq.pdf",
        cibersort="results/immune/cibersort.pdf",
        aggregate="results/immune/aggregate.pdf",
        rds="results/immune/immunedeconv.rds"
    log:
        "logs/immunedeconv.log",
    params:
        cibersort_path=config["immune"]["cibersort_path"],
        mcp_path=config["immune"]["mcp_path"],
        conda=config['env']['conda_shell'],
        species=config['ref']['species'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    shell:
        """
        source {params.conda} && conda activate {params.env};
        
        Rscript scripts/immunedeconv.R \
        --log {log} \
        --cibersort {output.cibersort_path} \
        --mcp {output.mcp_path} \
        --species {params.species} \
        --input {input} \
        --output {output.rds}
        """
