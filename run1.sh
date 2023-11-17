#!/bin/bash

# Add bioconda to conda channels
conda config --add channels bioconda

# Define your genome file
genome_file="genome.file.fasta"

# Define your chromosomes
reference_chromosomes="Chr01.2 Chr02.1 Chr03.1 ..."
alternative_chromosomes="Chr01.1 Chr02.2 Chr03.2 ..."

# Create reference and alternative fasta files
samtools faidx $genome_file $reference_chromosomes > Genome.reference.fasta
samtools faidx $genome_file $alternative_chromosomes > Genome.alternative.fasta

# Create and activate conda environment
conda create --name repeats
conda activate repeats

# Install RepeatModeler and RepeatMasker
conda install -c bioconda repeatmodeler
conda install -c bioconda repeatmasker

# Run RepeatModeler
BuildDatabase -name species_name $genome_file
nohup RepeatModeler -database species_name -pa 25 -LTRStruct &

# Run RepeatMasker
nohup RepeatMasker -gff -a -pa 20 -u -lib repeatmodeler_output-families.fa $genome_file &

# Mask repetitive regions
bedtools maskfasta -fi Genome.reference.fasta -bed out_put_from_RepeatMasker.fasta.out.gff -soft -fo Genome.reference.softmask.fasta
