# Functional Annotation Pipeline

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)


## Prerequisites

- **Python**: Version 3.6 or higher
- **Snakemake**: Workflow management system
- **Singularity**: Containerization platform

### Singularity Images

The necessary `.sif` and `.img` Singularity images are hosted on [Sylabs Cloud](https://cloud.sylabs.io/library/danifilho/functional_annotation_images). You can download them using the following commands:

```bash
singularity pull --arch amd64 library://danifilho/structural_annotation/validate_gff:sha256.e438afcf3e0e3c3e892adeb8bc85fa04d93c6616031780cc80e4e6177f99d4cf -F images/validate_gff.sif
singularity pull --arch amd64 library://danifilho/structural_annotation/rename_gff:sha256.98c0b56f2adef78a924496d9611a1ed56399d1a4fa0f213e077c4d6d65f5d9dc -F images/gffread.sif
singularity pull --arch amd64 library://danifilho/structural_annotation/custom_python:sha256.1a24456dca10a0b791d929fec50a7db0d768a9e340b31f0134d2053a2e160cd1 -F images/python.sif
singularity pull --arch amd64 library://danifilho/structural_annotation/maker:sha256.60d7f06b2d2cc97fd764b5002d9873dd14b1eacbf6f58645515e47416ee5828b -F images/maker.sif
singularity pull docker://greensii_isoseq3 
singularity pull --arch amd64 library://danifilho/structural_annotation/bedtools:sha256.2a0840734f789467ee1a2357a55fd23db4c4eac8ee41adb39902900ed45ca07a -F images/bedtools.sif
singularity pull --arch amd64 library://danifilho/functional_annotation_images/samtools:sha256.756b3e649207b774365c7e35edcbe106b644e345baeb3f2aee77285e1a4799be -F images/samtools.sif
singularity pull docker://dfam_tetools -F images/dfam_tetools.sif
```

## Installation

Basically, to run the structural annotation step you need:
- One config.yaml file, with your variables and inputs
- One Snakefile
- One "inputs" folder, with the genome in .fasta format, the .faa protein file and SMRT Pacbio RNA reads (usually .fastq)
- One "scripts" folder, with the following scripts: maker_control_files, validate_gff.py, renameGff.py, and gff_filter.py
- One "images" folder, with all the above singularity images. 


4. **Ensure Singularity is installed**:

   - Install Singularity following the official [Singularity Installation Guide](https://sylabs.io/guides/3.0/user-guide/installation.html).

## Dependencies

### Software

- **Snakemake**: Workflow management system
- **Singularity**: Containerization platform

### Python Packages

- **biopython**

### Running the Pipeline

```bash
snakemake --cores <number_of_cores> --use-singularity
```
## License

(LICENSE)

## Contact
dasilvaf@msu.edu

For any questions or issues, please contact:

- **idk yet what to put here**
- **Email**: email@example.com
- **GitHub**: [username](https://github.com/username)
