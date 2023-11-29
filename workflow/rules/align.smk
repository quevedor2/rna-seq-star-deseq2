rule align_pe:
  input:
    fq1=get_fq1,
    fq2=get_fq2
  output:
    alignedcoord="results/star/pe/{sample}/Aligned.sortedByCoord.out.bam",
    alignedtranscriptome="results/star/pe/{sample}/Aligned.toTranscriptome.out.bam",
  threads: 16
  params:
    stargtf=config['star']['gtf'],
    starparams=config['star']['params'],
    stargenome=config["star"]["star-genome"],
    outprefix="results/star/pe/HPB-070-merged/",
  shell:
    """
    module load STAR/2.7.9a
    
    STAR \
    --quantMode GeneCounts TranscriptomeSAM \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterIntronMotifs RemoveNoncanonical \
    --chimSegmentMin 10 \
    --chimOutType SeparateSAMold \
    --outSAMunmapped Within \
    --sjdbGTFfile {params.stargtf} {params.starparams} \
    --runThreadN {threads} \
    --genomeDir {params.stargenome} \
    --readFilesIn {input.fq1} {input.fq2} \
    --readFilesCommand zcat \
    --outFileNamePrefix {params.outprefix} \
    --outStd Log
    """
    
    
# rule align_se:
#     input:
#         fq1=get_map_reads_input_R1,
#     output:
#         "results/star/se/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
#         "results/star/se/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
#         "results/star/se/{sample}-{unit}/ReadsPerGene.out.tab",
#     log:
#         "logs/star-se/{sample}-{unit}.log",
#     params:
#         index=config["star"]["star-genome"],
#         extra="--quantMode GeneCounts TranscriptomeSAM "
#         "--outSAMtype BAM SortedByCoordinate "
#         "--outFilterIntronMotifs RemoveNoncanonical "
#         "--chimSegmentMin 10 "
#         "--chimOutType SeparateSAMold "
#         "--outSAMunmapped Within "
#         "--sjdbGTFfile {} {}".format(
#             config["star"]["gtf"], config["star"]["params"]
#         ),
#     threads: 16
#     wrapper:
#         "v0.75.0/bio/star/align"

rule index_coord:
  input:
    "results/star/pe/{sample}/Aligned.sortedByCoord.out.bam",
  output:
    "results/star/pe/{sample}/Aligned.sortedByCoord.out.bam.bai",
  params:
  shell:
    """
    module load samtools/1.17
    
    samtools index ${input} ${output}
    """
