#!/bin/bash
#SBATCH -J rnaseq
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH -t 10:00:00
#SBATCH -o rnaseq_F_%j.out
#SBATCH -e rnaseq_F_%j.err

# write this script to stdout-file
cat $0

# load module
module load Nextflow/23.10.0

# set the project directory
cd /home/materna/lu2023-12-14/Students/Jakob/merged/F

# run the command
nextflow run nf-core/rnaseq \
    --input /home/materna/lu2023-12-14/Students/Jakob/merged/F/samplesheet_F.csv \
    --outdir /home/materna/lu2023-12-14/Students/Jakob/merged/F/nfcore \
    --fasta /home/materna/lu2023-12-14/Students/Jakob/merged/reference/merged/merged.fa \
    --gtf /home/materna/lu2023-12-14/Students/Jakob/merged/reference/merged/merged.gtf \
    --skipQC \
    --skip_biotype_qc \
    --skip_markduplicates \
    --skip_bigwig \
    --skip_stringtie \
    --skip_preseq \
    --skip_dupradar \
    --skip_qualimap \
    --skip_rseqc \
    --skip_deseq2_qc \
    --skip_multiqc \
    --max_cpus 48 \
    --max_memory 254GB \
    -profile singularity 
    