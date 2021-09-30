rule prepare_reference:
  input:
    reference_genome=config['ref_index']['genome'],
  output:
    seq="ref/reference.seq",
    grp="ref/reference.grp",
    ti="ref/reference.ti",
  params:
    extra="--gtf {}".format(config["star"]["gtf"]),
  log:
    "logs/rsem/prepare-reference.log",
  wrapper:
    "0.77.0/bio/rsem/prepare-reference"

rule calculate_expression:
  input:
    bam=get_star_transcriptome,
    reference="ref/reference.seq",
  output:
    genes_results="results/rsem/{sample}-{unit}.genes.results",
    isoforms_results="results/rsem/{sample}-{unit}.isoforms.results",
    tbam=temp("results/rsem/{sample}-{unit}.transcript.bam"),
    gbam=temp("results/rsem/{sample}-{unit}.genome.bam"),
  params:
    paired_end=lambda w: is_paired_end(w.sample),
    extra="-bam --estimate-rspd --output-genome-bam --time --forward-prob 0 --seed 42",
  log:
    "logs/rsem/calculate_expression/{sample}-{unit}.log",
  wrapper:
    "0.77.0/bio/rsem/calculate-expression"

rule rsem_generate_data_matrix:
  input:
    get_rsem_output_all_units,
  output:
    temp("results/counts/all.tmp"),
  params:
    extra="",
  log:
    "logs/rsem/generate_data_matrix.log",
  wrapper:
    "0.77.0/bio/rsem/generate-data-matrix"

rule format_data_matrix:
  input:
    "results/counts/all.tmp",
  output:
    "results/counts/all.tsv",
  shell:
    "sed 's/\"//g' {input} |  "
    "sed 's/results\/rsem\///g' | "
    "sed 's/-merged.genes.results//g' | "
    "sed '1 s/^/gene/' | "
    "sed -e 's/\.[0-9]*//g' > {output}"
