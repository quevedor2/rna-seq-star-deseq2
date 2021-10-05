import os
import sys
from snakemake.shell import shell

# the reference_name argument is inferred by stripping the .seq suffix from
# the output.seq value
# sys.argv[1] = snakemake.output.seq
output_directory = os.path.dirname(os.path.abspath(sys.argv[1]))
seq_file = os.path.basename(sys.argv[1])
if seq_file.endswith(".seq"):
    reference_name = os.path.join(output_directory, seq_file[:-4])
else:
    raise Exception("output.seq has an invalid file suffix (must be .seq)")

# sys.argv[2] = snakemake.params.get("extra", "")
extra = sys.argv[2]
# sys.argv[3] = snakemake.log_fmt_shell(stdout=True, stderr=True)
log = sys.argv[3]
shell(
    "rsem-prepare-reference --num-threads {snakemake.threads} {extra} "
    "{snakemake.input.reference_genome} {reference_name} > "
    "{log} 2>1"
)
