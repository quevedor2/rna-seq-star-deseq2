rule replace_rg:
  input:
    get_star_bam
  output:
    "results/star/rg/{sample}/Aligned.rg.bam",
  log:
    "logs/picard/replace_rg/{sample}.log"
  params:
    "RGLB=lib1 RGPL=illumina RGPU={sample} RGSM={sample}"
  resources:
    mem_mb=1024
  wrapper:
    "v0.75.0/bio/picard/addorreplacereadgroups"

rule index_rg:
  input:
    get_rg_bam,
  output:
    "results/star/rg/{sample}/Aligned.rg.bam.bai",
  log:
    "logs/samtools/index/{sample}.log"
  wrapper:
    "v0.75.0/bio/samtools/index"

rule collect_allelic_counts:
  input:
    bam=get_rg_bam,
    bai=get_rg_bai,
  output:
    "results/genotyping/{sample}.allelicCounts.tsv",
  conda:
    "../envs/gatk.yaml",
  log:
    "logs/gatk/collectalleliccounts/{sample}.log"
  params:
    ref=config['ref_index']['genome'],
    target=config['genotyping']['target'],
  shell:
    "touch {output}; "
    "echo {input.bam} > {log}; "
    "echo {input.bai} >> {log}; "
    "gatk --java-options '-Xmx59G  -XX:ParallelGCThreads=1' CollectAllelicCounts "
    "-I {input.bam} "
    "-R {params.ref} "
    "-L {params.target} "
    "-O {output} > {log} 2>&1"

rule categorizeAD_gatk:
  input:
    "results/genotyping/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/genotyping/{sample}_out.tmp"),
    simple=temp("results/genotyping/{sample}_out.tsv"),
  params:
    ref=2,
    alt=3,
    mincov=10,
  conda:
    "../envs/perl.yaml",
  shell:
    "grep -v '^@' {input} | tail -n +2 > {output.intermediate}; "
    "perl scripts/allelic_count_helper.pl categorize "
    "{output.intermediate} {params.ref} {params.alt} {params.mincov} > {output.simple}"


rule aggregate_AD:
  input:
    expand(
        "results/genotyping/{unit.sample_name}-{unit.unit_name}_out.tsv",
        unit=units.itertuples(),
    ),
  output:
    "results/zygosity/AD/aggregate.csv",
  params:
    ""
  conda:
    "../envs/perl.yaml",
  shell:
    "paste -d',' {input} > {output}; "

rule get_AD_lines:
  input:
    "results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_lines.txt",
  params:
    min_n=2,
  conda:
    "../envs/perl.yaml",
  shell:
    "perl scripts/allelic_count_helper.pl setlines "
    "{input} {params.min_n} > {output}; "

rule filt_AD_raw:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
    aggregate_raw="results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_filt.csv",
  params:
    samples=expand(
        "{unit.sample_name}-{unit.unit_name}",
        unit=units.itertuples(),
    ),
  conda:
    "../envs/perl.yaml",
  shell:
    "echo {params.samples} | sed 's/\s/,/g' > {output}; "
    "perl scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {input.aggregate_raw} >> {output}; "

rule filt_AD_dbsnp:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
  output:
    "results/zygosity/AD/aggregate_pos.txt",
  params:
    target=config['genotyping']['target'],
  conda:
    "../envs/perl.yaml",
  shell:
    "perl scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {params.target} > {output}; "

rule run_wp_zygosity:
  input:
    het_pos='results/zygosity/AD/aggregate_pos.txt',
    het_cnt='results/zygosity/AD/aggregate_filt.csv',
  output:
    pdfs=report(expand("results/zygosity/wadingpool/hmmfit_{sample}.pdf", sample=samples.index),
                caption="../report/wp_zygosity.rst", category="CNV"),
    tbls=expand("results/zygosity/wadingpool/hmmfit_{sample}.tsv", sample=samples.index),
    rdas=expand("results/zygosity/wadingpool/hmmfit_{sample}.rda", sample=samples.index),
  params:
    outdir='results/zygosity/wadingpool',
    cndir='results/cnv/ichorcna',
    genome=config['ref']['build'],
    maxstate=300,
  conda:
    "../envs/r.yaml",
  shell:
    "wp_zygosity.R "
    "--hetposf {input.het_pos} "
    "--hetf {input.het_cnt} "
    "--cnpath {params.cndir} "
    "--outdir {params.outdir} "
    "--genome {params.genome} "
    "--maxstate {params.maxstate} "

rule run_wp_identity_autosome:
  input:
    het_cnt='results/zygosity/AD/aggregate_filt.csv',
  output:
    pdf=report("results/sampleid/autosome/similarity_autosome.pdf",
                caption="../report/autoid.rst", category="genotypeID"),
    tbl=report("results/sampleid/autosome/similarity_autosome.tsv",
                caption="../report/autoid.rst", category="genotypeID"),
  params:
    outdir='results/sampleid/autosome',
    midpoint=0.03,
  conda:
    "../envs/r.yaml",
  shell:
    "wp_identity.R "
    "--sample {input.het_cnt} "
    "--midpoint {params.midpoint} "
    "--mode autosome "
    "--outdir {params.outdir} "
