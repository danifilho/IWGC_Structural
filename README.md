# Functional Annotation Pipeline

## Table of Contents

- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Configuration](#configuration)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)


## Prerequisites

- **Python**: Version 3.6 or higher
- **Snakemake**: Workflow management system
- **SingularityCE**: Containerization platform (CE version works better with docker/singularity)

### Singularity Images

The necessary `.sif` Singularity images are hosted on [Sylabs Cloud](https://cloud.sylabs.io/library/danifilho/functional_annotation_images). You can download them using the following commands:

```bash
# Structural annotation images (Sylabs Library - unsigned => use -U)
singularity pull -U -F validate_gff.sif  library://danifilho/structural_annotation/validate_gff:latest
singularity pull -U -F rename_gff.sif    library://danifilho/structural_annotation/rename_gff:latest
singularity pull -U -F python.sif        library://danifilho/structural_annotation/custom_python:latest
singularity pull -U -F maker.sif         library://danifilho/structural_annotation/maker:latest
singularity pull -U -F bedtools.sif      library://danifilho/structural_annotation/bedtools:latest

# IsoSeq3 (Docker registry)
singularity pull -F greensii_isoseq3.sif docker://greensii/isoseq3

# Functional annotation images (Sylabs Library - unsigned => use -U)
singularity pull -U -F samtools.sif      library://danifilho/functional_annotation_images/samtools:latest

# TE tools (Docker registry)
singularity pull -F dfam_tetools.sif     docker://dfam/tetools
```

## Setup

Basically, to run the structural annotation step you need:
- One config.yaml file, with your variables and inputs
- One Snakefile
- One "inputs" folder, with the genome in .fasta format, the .faa protein file and SMRT Pacbio RNA reads (usually .fastq)
- One "scripts" folder, with the following scripts: maker_control_files, validate_gff.py, renameGff.py, and gff_filter.py
- One "images" folder, with all the above singularity images.

## Configuration

Create a `config.yaml` file in the root directory with the following content (example), or copy from this github:

```yaml
fasta_file: chromosome_1.fasta
isoseq3_reads: SRR28834784.fastq
protein_file: protein.faa
volume_name: /mnt/ufs18/rs-028/Weed_lab/software/testing/IWGC_annotation/structural
species_abbreviation: BetVu

#REPEATING PARAMS
repeatmodeler_threads: 90
repeatmasker_threads: 90

#ISOSEQ3 PARAMS
isoseq3_min_coverage: 0.9
isoseq3_min_identity: 0.95

#PBMM2 PARAMS
pbmm2_threads: 40
```

## Dependencies

### Software

- **Snakemake**: Workflow management system
- **Singularity**: Containerization platform

### Python Packages

- **biopython**

## Usage

```bash
snakemake --cores <number_of_cores> --use-singularity
```
## License

(LICENSE)

## Contact

- **Email**: dasilvaf@msu.edu
- **GitHub**: [danifilho](https://github.com/danifilho)
