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
        renv=config['env']['renv_path'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    shell:
        """
        module load R/4.2.1
        
        Rscript scripts/immunedeconv.R \
        --log {log} \
        --renv {params.renv} \
        --cibersort {params.cibersort_path} \
        --mcp {params.mcp_path} \
        --species {params.species} \
        --input {input} \
        --output {output.rds}
        """
