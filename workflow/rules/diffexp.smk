rule count_matrix:
    input:
        get_star_output_all_units,
    output:
        "results/counts/all_old.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="results/counts/all.tsv",
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        """
        source {params.conda} && conda activate {params.env};
        
        ../scripts/deseq2-init.R
        """


rule pca:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/pca.pdf", "../report/pca.rst"),
    params:
        pca_labels=config["pca"]["labels"],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    log:
        "logs/pca.log",
    script:
        """
        source {params.conda} && conda activate {params.env};
        
        ../scripts/plot-pca.R
        """


rule deseq2:
    input:
        "results/deseq2/all.rds",
    output:
        table=report(
            "results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"
        ),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.pdf", "../report/ma.rst"),
    params:
        contrast=get_contrast,
        species=config['ref']['species'],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    script:
        """
        source {params.conda} && conda activate {params.env};
        
        ../scripts/deseq2.R
        """

rule gsea:
    input:
        "results/diffexp/{contrast}.diffexp.tsv",
    output:
        ora_tbl=report("results/gsea/{contrast}.ora.tsv", "../report/ora-go.rst"),
        gsea_tbl=report("results/gsea/{contrast}.gsea.tsv", "../report/gsea-kegg.rst"),
        gsea_pdf="results/gsea/{contrast}.gsea.pdf",
    params:
        genome=config['ref']['species'],
        contrast=get_contrast,
        minbase=config["gsea"]["min_base_mean"],
        maxp=config["gsea"]["max_p"],
        conda=config['env']['conda_shell'],
        env=directory(config['env']['r41']),
#    conda:
#        "../envs/deseq2.yaml"
    log:
        "logs/gsea/{contrast}.log",
    threads: get_deseq2_threads
    script:
        """
        source {params.conda} && conda activate {params.env};
        
        ../scripts/gsea.R
        """
