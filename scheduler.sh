#! /bin/bash
#SBATCH -J scheduler
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p all
#SBATCH --mem=2gb

cd /cluster/home/quever/workflows/rna-seq-star-deseq2
condaprefix='/cluster/home/quever/workflows/rna-seq-star-deseq2/.snakemake/conda'

snakemake \
--jobs 6 \
--profile slurm \
--cluster-config slurm/cluster.json \
--conda-frontend conda \
--use-conda \
--use-singularity \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:///cluster/home/quever/downloads/snakemake-wrappers/' \
--rerun-incomplete
