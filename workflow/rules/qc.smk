## RSEQC
rule rseqc_gtf2bed:
    input:
        config["star"]["gtf"],
    output:
        bed="results/qc/rseqc/annotation.bed",
    log:
        "logs/rseqc_gtf2bed.log",
    shell:
        """
        module load  ucsctools/378
        
        cat {input} |\
          gtfToGenePred /dev/stdin /dev/stdout |\
          genePredToBed /dev/stdin /dev/stdout \
          > {output.bed}
        """


rule fastqcWIP:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/fastqc/{sample}.merged_fastqc.html",
    priority: 1
    log:
        "logs/fastqc/{sample}.log",
    params:
        outdir="results/qc/fastqc/",
    shell:
        """
        module load fastqc/0.11.5
        
        outputdir=$(dirname {output})
        
        fastqc \
        -f bam \
        -o $outputdir \
        {input.bam}
        """

rule rseqc_junction_annotation:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            junction_annotation.py {params.extra} \
            -i {input.bam} \
            -r {input.bed} \
            -o {params.prefix} \
            > {log[0]} 2>&1
        """


rule rseqc_junction_saturation:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            junction_saturation.py {params.extra} \
            -i {input.bam} \
            -r {input.bed} \
            -o {params.prefix}
            > {log} 2>&1
        """


rule rseqc_stat:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
    output:
        "results/qc/rseqc/{sample}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}.log",
    params:
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            bam_stat.py \
            -i {input.bam} > {output} 2> {log}
        """


rule rseqc_infer:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}.log",
    params:
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}
        """


rule rseqc_innerdis:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            inner_distance.py -r {input.bed} \
            -i {input.bam} -o {params.prefix} > {log} 2>&1
        """


rule rseqc_readdis:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}.log",
    params:
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            read_distribution.py \
            -r {input.bed} \
            -i {input.bam} > {output} 2> {log}
        """


rule rseqc_readdup:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
    output:
        "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            read_duplication.py \
            -i {input.bam} \
            -o {params.prefix} > {log} 2>&1
        """


rule rseqc_readgc:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
    output:
        "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            read_GC.py \
            -i {input.bam} \
            -o {params.prefix} > {log} 2>&1
        """

rule rseqc_readnvc:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
    output:
        "results/qc/rseqc/{sample}.read_nvc.NVC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readnvc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".NVC_plot.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            read_NVC.py \
            -i {input.bam} \
            -o {params.prefix} > {log} 2>&1
        """


rule rseqc_rpkmsaturation:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_rpkmsaturation.saturation.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_rpkmsaturation/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".saturation.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            RPKM_saturation.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.prefix} > {log} 2>&1
        """

rule rseqc_genebodycoverage:
    input:
        bam=lambda wc: get_star_output(wc, fi='coord'),
        bai=lambda wc: get_star_output(wc, fi='coord', bai=True),
    output:
        "results/qc/rseqc/{sample}_genebodycoverage.geneBodyCoverage.curves.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_genebodycoverage/{sample}.log",
    params:
        housekeeping=config["rseqc"]["housekeeping_genes"],
        prefix=lambda w, output: output[0].replace(".geneBodyCoverage.curves.pdf", ""),
        rseqc=config['env']['rseqc']
    shell:
        """
        module load apptainer/1.0.2
        
        apptainer exec {params.rseqc} \
            geneBody_coverage.py \
            -r {params.housekeeping} \
            -i {input.bam} \
            -o {params.prefix} > {log} 2>&1
        """

rule multiqc:
    input:
        lambda wc: get_star_output_all_units(wc, fi="bam"),
        expand(
            "results/qc/rseqc/{unit.sample_name}.junctionanno.junction.bed",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.infer_experiment.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.stats.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.readdistribution.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "logs/rseqc/rseqc_junction_annotation/{unit.sample_name}.log",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}.read_nvc.NVC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "results/qc/rseqc/{unit.sample_name}_rpkmsaturation.saturation.pdf",
            unit=units.itertuples(),
        ),
        expand(
           "results/qc/rseqc/{unit.sample_name}_genebodycoverage.geneBodyCoverage.geneBodyCoverage.curves.pdf",
           unit=units.itertuples(),
        ),
        expand(
           "results/qc/fastqc/{unit.sample_name}.merged_fastqc.html",
           unit=units.itertuples(),
        ),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    params:
        multiqc=config['env']['multiqc'],
        rseqcdir='results/qc/rseqc/',
        fastqcdir='results/qc/fastqc/',
    shell:
        """
        module load apptainer/1.0.2
        
        outputdir=$(dirname {output})
        outputfile=$(basename {output})
        
        apptainer exec {params.multiqc} \
            multiqc \
            --force \
            -o $(dirname {output}) \
            -n $(basename {output}) \
            {params.rseqcdir} {params.fastqcdir}
        """
