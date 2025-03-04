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
- **Singularity**: Containerization platform

### Singularity Images

The necessary `.sif` Singularity images are hosted on [Sylabs Cloud](https://cloud.sylabs.io/library/danifilho/functional_annotation_images). You can download them using the following commands:

```bash
singularity pull validate_gff.sif library://danifilho/structural_annotation/validate_gff:sha256.e438afcf3e0e3c3e892adeb8bc85fa04d93c6616031780cc80e4e6177f99d4cf
singularity pull gffread.sif library://danifilho/structural_annotation/rename_gff:sha256.98c0b56f2adef78a924496d9611a1ed56399d1a4fa0f213e077c4d6d65f5d9dc 
singularity pull python.sif library://danifilho/structural_annotation/custom_python:sha256.1a24456dca10a0b791d929fec50a7db0d768a9e340b31f0134d2053a2e160cd1
singularity pull maker.sif library://danifilho/structural_annotation/maker:sha256.60d7f06b2d2cc97fd764b5002d9873dd14b1eacbf6f58645515e47416ee5828b
singularity pull greensii_isoseq3.sif docker://greensii/isoseq3
singularity pull bedtools.sif library://danifilho/structural_annotation/bedtools:sha256.2a0840734f789467ee1a2357a55fd23db4c4eac8ee41adb39902900ed45ca07a
singularity pull samtools.sif library://danifilho/functional_annotation_images/samtools:sha256.756b3e649207b774365c7e35edcbe106b644e345baeb3f2aee77285e1a4799be
singularity pull dfma_tetools.sif docker://dfam/tetools 
```

## Setup

Basically, to run the structural annotation step you need:
- One config.yaml file, with your variables and inputs
- One Snakefile
- One "inputs" folder, with the genome in .fasta format, the .faa protein file and SMRT Pacbio RNA reads (usually .fastq)
- One "scripts" folder, with the following scripts: maker_control_files, validate_gff.py, renameGff.py, and gff_filter.py
- One "images" folder, with all the above singularity images.

## Configuration

Create a `config.yaml` file in the root directory with the following content:

```yaml
species_name: "your_species_name"
volume_name: "/absolute/path/to/your/working/directory"
fasta_file: "genome.fasta"
final_filtering_gff: "annotations.gff"
multiloc_script: "multiloc_script.py"
mmseqs_databases:
  - name: "database1"
    path: "/absolute/path/to/databases/database1"
  - name: "ncbi"
    path: "/absolute/path/to/ncbi_data_file.prot"
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
