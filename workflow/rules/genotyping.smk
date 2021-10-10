rule replace_rg_pe:
  input:
    get_star_bam,
  output:
    "results/star/pe/{sample}-{unit}/Aligned.rg.bam",
  log:
    "logs/picard/replace_rg_pe/{sample}-{unit}.log"
  params:
    "RGLB=lib1 RGPL=illumina RGPU={sample}-{unit} RGSM={sample}-{unit}"
  resources:
    mem_mb=1024
  wrapper:
    "v0.75.0/bio/picard/addorreplacereadgroups"

rule replace_rg_se:
  input:
    get_star_bam,
  output:
    "results/star/se/{sample}-{unit}/Aligned.rg.bam",
  log:
    "logs/picard/replace_rg_se/{sample}-{unit}.log"
  params:
    "RGLB=lib1 RGPL=illumina RGPU={sample}-{unit} RGSM={sample}-{unit}"
  resources:
    mem_mb=1024
  wrapper:
    "v0.75.0/bio/picard/addorreplacereadgroups"

rule collect_allelic_counts:
  input:
    bam=get_rg_bam,
    bai=get_rg_bai,
  output:
    "results/genotyping/{sample}-{unit}.allelicCounts.tsv",
  conda:
    "../envs/gatk.yaml",
  log:
    "logs/gatk/collectalleliccounts/{sample}-{unit}.log"
  params:
    ref=config['ref_index']['genome'],
    target=config['genotyping']['target'],
  shell:
    "touch {output}; "
    "echo {input.bam} > {log}; "
    "echo {input.bai} >> {log}; "
    "gatk --java-options '-Xmx20G  -XX:ParallelGCThreads=4' CollectAllelicCounts "
    "-I {input.bam} "
    "-R {params.ref} "
    "-L {params.target} "
    "-O {output} > {log} 2>&1"

rule categorizeAD_gatk:
  input:
    "results/zygosity/counts/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/zygosity/counts/{sample}_out.tmp"),
    simple=temp("results/zygosity/counts/{sample}_out.tsv"),
  params:
    ref=2,
    alt=3,
  conda:
    "../envs/perl.yaml",
  shell:
    "grep -v '^@' {input} | tail -n +2 > {output.intermediate}; "
    "perl workflow/scripts/allelic_count_helper.pl categorize "
    "{output.intermediate} {params.ref} {params.alt} > {output.simple}"


rule aggregate_AD:
  input:
    expand("results/zygosity/counts/{sample}_out.tsv", sample=samples.index),
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
    "perl workflow/scripts/allelic_count_helper.pl setlines "
    "{input} {params.min_n} > {output}; "

rule filt_AD_raw:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
    aggregate_raw="results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_filt.csv",
  params:
    samples=expand("{sample}", sample=samples.index),
  conda:
    "../envs/perl.yaml",
  shell:
    "echo {params.samples} | sed 's/\s/,/g' > {output}; "
    "perl workflow/scripts/allelic_count_helper.pl getlines "
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
    "perl workflow/scripts/allelic_count_helper.pl getlines "
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
