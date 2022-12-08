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
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    script:
        """
        source {params.conda} && conda activate {params.env};
        
        ../scripts/immunedeconv.R
        """
