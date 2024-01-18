import glob
import pandas as pd
from snakemake.utils import validate
#from snakemake.remote import FTP
#ftp = FTP.RemoteProvider()
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)


def get_final_output():
    final_output = expand(
        "results/gsea/{contrast}.gsea.pdf",
        contrast=config["diffexp"]["contrasts"],
    )
    return final_output


validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")



def get_fq1(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_paired_end(wildcards.sample):
        u = units.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
    else:
        u = units.loc[ (wildcards.sample), ["fq1"] ].dropna()
    return [ f"{u.fq1}" ]

def get_fq2(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_paired_end(wildcards.sample):
        u = units.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
        fq2 = [ f"{u.fq2}" ]
    else:
        fq2 = ""
    return fq2

def get_rsem_output_all_units(wildcards):
    res = []
    for unit in units.itertuples():
        res.append(
            "results/rsem/{}.genes.results".format(
                unit.sample_name
            )
        )
    return res

def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units.isnull()["fq2"]
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired
    
def get_star_output_all_units(wildcards, fi="coord", use='all'):
    if fi == "coord":
        outfile = "Aligned.sortedByCoord.out.bam"
    elif fi == 'transcriptome':
        outfile = "Aligned.toTranscriptome.out.bam"
    else:
	outfile = "ReadsPerGene.out.tab"
    res = []
    if use == 'all':
        for unit in units.itertuples():
            if is_paired_end(unit.sample_name):
                lib = "pe"
            else:
                lib = "se"
            res.append(
                "results/star/{}/{}/{}".format(
                    lib, unit.sample_name, outfile
                )
            )
    else:
        if is_paired_end(wildcards.sample):
            lib = 'pe'
        else:
            lib = 'se'
        res.append("results/star/{}/{}/{}".format(
                lib, wildcards.sample, outfile
                )
        )
    return res

def get_star_output(wildcards, fi="coord", bai=False):
    if fi == "coord":
        outfile = "Aligned.sortedByCoord.out.bam"
    elif fi == 'transcriptome':
        outfile = "Aligned.toTranscriptome.out.bam"
    else:
        outfile = "ReadsPerGene.out.tab"
    if bai:
        outfile = outfile + ".bai"
        
    if is_paired_end(wildcards.sample):
        lib = 'pe'
    else:
        lib = 'se'
    
    res = []
    res.append("results/star/{}/{}/{}".format(
        lib, wildcards.sample, outfile
        )
    )
    return res
    
