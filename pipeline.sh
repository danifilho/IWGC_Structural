
#!/bin/bash
########## SBATCH Lines for Resource Request ##########


# Load the configuration file

module load GCC/8.3.0
module load OpenMPI/3.1.4
module load RepeatModeler
module load RepeatMasker
module load BEDTools

# build variable names from species name
FASTA_FILE=$1
DATABASE_NAME=$2
NUM_JOBS_REPEAT=$3
RMODELER_OUTPUT=$4
RMASKER_OUTPUT=$5
BEDTOOLS_OUTPUT=$6
ISOSEQ3_READS=$7
PBMM2_OUTPUT=$8
ISOSEQ3_OUTPUT=$9
#--------------------------- REPEAT MODELER

# Use the database name and fasta file from the configuration file
BuildDatabase -name $DATABASE_NAME $FASTA_FILE

# Use the database name and number of jobs (threads) from the configuration file
RepeatModeler -database $DATABASE_NAME -pa $NUM_JOBS_REPEAT -LTRStruct
#---------------------------- REPEAT MASKER
## Use the output files as input for the next command
RepeatMasker -gff -a -pa 20 -u $RMODELER_OUTPUT $FASTA_FILE
#----------------------------- BEDTOOLS

bedtools maskfasta -fi $FASTA_FILE -bed $RMASKER_OUTPUT -soft -fo $BEDTOOLS_OUTPUT

#----------------------------- MAPPING
#cd /mnt/gs21/scratch/cuttilua/BASF/Col_PEK1.5_assembly_and_annotation/Mapping

#source /mnt/home/cuttilua/anaconda3/bin/activate

#conda activate basf

$ module --ignore_cache load "isoseq3"
module loadpbmm2

pbmm2 align --preset ISOSEQ -O6,24 -B4 --sort $ISOSEQ3_READS $BEDTOOLS_OUTPUT $PBMM2_OUTPUT -j 36 #WE SHOULD CHANGE IT LATER

isoseq3 collapse --do-not-collapse-extra-5exons --min-aln-coverage 0.9 --min-aln-identity 0.95 $PBMM2_OUTPUT $ISOSEQ3_OUTPUT

