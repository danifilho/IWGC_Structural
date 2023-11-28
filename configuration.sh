#!/bin/bash

#SBATCH --time=36:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=36           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=300G           # memory required per node. can also make mem per cpu with different --mem-per-cpu=2G
#SBATCH --job-name EDTA_GR     # you can give your job a name for easier identification (same as -J)

# Configuration for RepeatModeler script

#Setting species name
SPECIES_NAME="Arabidopsis_thaliana"

# Path to RepeatModeler
REPEATMODELER_PATH="/opt/software/RepeatModeler/2.0.2a-foss-2019b/RepeatModeler"

# Name of the database
DATABASE_NAME="{$SPECIES_NAME}_DB"

# Name of the fasta file
FASTA_FILE="Arabidopsis_thaliana.fasta"

# Number of parallel jobs
NUM_JOBS_REPEAT=36

# Find the repeatmodeler output file
RMODELER_OUTPUT=$SPECIES_NAME-families.fa

# Using repeatmasker output
RMASKER_OUTPUT=$SPECIES_NAME.fasta.out.gff

#Using bedtools output
BEDTOOLS_OUTPUT=$SPECIES_NAME.softmask.fasta

#------------------------------- mapping
ISOSEQ3_READS="AraTh_Leaves_SRR14584396.fastq"

PBMM2_OUTPUT=$SPECIES_NAME-aligned.bam
ISOSEQ3_OUTPUT=$SPECIES_NAME.collapsed.gff

#modify by 20/11/2023
./folders.sh 
./pipeline.sh $FASTA_FILE $DATABASE_NAME $NUM_JOBS_REPEAT $RMODELER_OUTPUT $RMASKER_OUTPUT $BEDTOOLS_OUTPUT $ISOSEQ3_READS $PBMM2_OUTPUT $ISOSEQ3_OUTPUT

./pipeline.sh 
