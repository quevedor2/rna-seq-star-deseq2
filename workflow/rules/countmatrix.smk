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
    bam="results/star/se/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    reference="ref/reference.seq",
  output:
    genes_results="results/rsem/{sample}-{unit}.genes.results",
    isoforms_results="results/rsem/{sample}-{unit}.isoforms.results",
  params:
    paired_end=lambda w: is_paired_end,
    extra="--seed 42",
  log:
    "logs/rsem/calculate_expression/{sample}-{unit}.log",
  wrapper:
    "0.77.0/bio/rsem/calculate-expression"

rule rsem_generate_data_matrix:
  input:
    get_rsem_output_all_units,
  output:
    "results/counts/all.tsv",
  params:
    extra="",
  log:
    "logs/rsem/generate_data_matrix.log",
  wrapper:
    "0.77.0/bio/rsem/generate-data-matrix"
