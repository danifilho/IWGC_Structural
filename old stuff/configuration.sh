#!/bin/bash

#====================== this configuration is for the HPCC line
#SBATCH --time=36:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=36           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=300G           # memory required per node. can also make mem per cpu with different --mem-per-cpu=2G
#SBATCH --job-name EDTA_GR     # you can give your job a name for easier identification (same as -J)

#Setting species name (It will be our reference to next outputs)
SPECIES_NAME="Arabidopsis_thaliana.Col-PEK1.5" #input

# Path to RepeatModeler
REPEATMODELER_PATH="/opt/software/RepeatModeler/2.0.2a-foss-2019b/RepeatModeler"

# Name of the database (necessary for repeatmodeler)
DATABASE_NAME=$SPECIES_NAME

# Name of the fasta file
FASTA_FILE="Arabidopsis_thaliana.Col-PEK1.5.fasta" #input

# Number of parallel jobs in REPEAT process (we can change it later)
NUM_JOBS_REPEAT=36

# Find the repeatmodeler output file
RMODELER_OUTPUT=$DATABASE_NAME-families.fa

# Using repeatmasker output
RMASKER_OUTPUT=$SPECIES_NAME.fasta.out.gff

#Using bedtools output
BEDTOOLS_OUTPUT=$SPECIES_NAME.softmask.fasta

#------------------------------- MAPPING ARGUMENTS
#isoseq3 input
ISOSEQ3_READS="AraTh_Leaves_SRR14584396.fastq" #input

# Using PBMM2 output
PBMM2_OUTPUT=$SPECIES_NAME-aligned.bam
# Using ISOSEQ3 output
ISOSEQ3_OUTPUT=$SPECIES_NAME.collapsed.gff

#-------------------------------- MAKER ARGUMENTS
# Input file for Maker
PROTEIN_FILE="Canola.PROTEINS.faa" #input

# Open the environment (to run in HPCC)

cd /mnt/gs21/scratch/cuttilua/BASF/Col_PEK1.5_assembly_and_annotation/Maker
source /mnt/home/cuttilua/anaconda3/bin/activate
conda activate basf

./pipeline.sh $FASTA_FILE $DATABASE_NAME $NUM_JOBS_REPEAT $RMODELER_OUTPUT $RMASKER_OUTPUT $BEDTOOLS_OUTPUT $ISOSEQ3_READS $PBMM2_OUTPUT $ISOSEQ3_OUTPUT $PROTEIN_FILE
