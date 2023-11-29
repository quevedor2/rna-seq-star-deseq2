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
    u = units.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq1}" ]

def get_fq2(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = units.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq2}" ]
