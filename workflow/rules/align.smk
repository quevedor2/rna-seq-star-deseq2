rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
    output:
        "results/star/pe/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/pe/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "results/star/pe/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
    params:
        index=config["star"]["star-genome"],
        extra="--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--chimSegmentMin 10 "
        "--chimOutType SeparateSAMold "
        "--outSAMunmapped Within "
        "--sjdbGTFfile {} {}".format(
            config["star"]["gtf"], config["star"]["params"]
        ),
    threads: 24
    wrapper:
        "v0.75.0/bio/star/align"

rule align_se:
    input:
        fq1=get_map_reads_input_R1,
    output:
        "results/star/se/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/se/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "results/star/se/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-se/{sample}-{unit}.log",
    params:
        index=config["star"]["star-genome"],
        extra="--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--chimSegmentMin 10 "
        "--chimOutType SeparateSAMold "
        "--outSAMunmapped Within "
        "--sjdbGTFfile {} {}".format(
            config["star"]["gtf"], config["star"]["params"]
        ),
    threads: 24
    wrapper:
        "v0.75.0/bio/star/align"

rule index_coord:
  input:
    get_star_bam,
  output:
    "results/star/{ends}/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
  log:
    "logs/samtools/index/{sample}-{unit}.{ends}.sortedByCoord.log"
  wrapper:
    "v0.75.0/bio/samtools/index"
