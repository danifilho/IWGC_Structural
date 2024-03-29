# Snakemake file: Snakefile
configfile: "config.yaml"

# Create a list of chromosomes from the fasta file
from Bio import SeqIO
chr = []
for record in SeqIO.parse(f"{config['volume_name']}/inputs/{config['fasta_file']}", "fasta"):
    chr.append(record.id)

# Take the fasta file name a save it as species_name, removing the fasta extension
if '.' in config['fasta_file']:
    species_name = config['fasta_file'][:config['fasta_file'].rindex('.')]  # Remove the extension
else:
    species_name = config['fasta_file']  # If no extension found, use the original file name

# Define the final output expected and the end of the pipeline
rule all:
    input:
        expand(f"filtering_outputs/{config['species_abbreviation']}.v2.gff")

# Create the database and identify repeated DNA elements in the genome
rule repeatmodeler:
    input:
        fasta_file=f"inputs/{config['fasta_file']}"
    output:
        repeatmodeler_output=f"repeat_modeler_outputs/{species_name}-families.fa"
    params:
        volume_name=config["volume_name"],
        num_jobs_repeat=config["num_jobs_repeat"],
        database_name=f"{species_name}"
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest BuildDatabase -name /data/repeat_modeler_outputs/{params.database_name} /data/{input.fasta_file} && "
        "docker run -v \"{params.volume_name}/repeat_modeler_outputs\":/data dfam/tetools:latest RepeatModeler -database /data/{params.database_name} -threads {params.num_jobs_repeat} -LTRStruct"

# Mask repetitive sequences
rule repeatmasker:
    input:
        fasta_file=f"inputs/{config['fasta_file']}",
        repeatmodeler_output=rules.repeatmodeler.output
    output:
        repeatmasker_output=f"repeat_masker_outputs/{species_name}.fasta.out.gff"
    params:
        volume_name=f"{config['volume_name']}",
        repeatmasker_threads=config["repeatmasker_threads"]
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatMasker -gff -a -pa {params.repeatmasker_threads} -u /data/repeat_modeler_output/{input.repeatmodeler_output} /data/{input.fasta_file} -dir /data/repeat_masker_outputs"

# Perform softmasking in repetitive sequences
rule bedtools:
    input:
        fasta_file=f"inputs/{config['fasta_file']}",
        repeatmasker_output=rules.repeatmasker.output
    output:
        bedtools_output=f"bedtools_outputs/{species_name}.softmask.fasta"
    params:
        volume_name=config["volume_name"]
    shell:
        "docker run -v \"{params.volume_name}\":/data pegi3s/bedtools:latest bedtools maskfasta -fi /data/{input.fasta_file} -bed /data/{input.repeatmasker_output} -soft -fo /data/{output.bedtools_output}"

# Align the isose reads with the masked genome
rule pbmm2:
    input:
        isoseq3_reads=f"inputs/{config['isoseq3_reads']}",
        bedtools_output=rules.bedtools.output
    output:
        pbmm2_output=f"pbmm2_outputs/{species_name}-aligned.bam"
    params:
        volume_name=config["volume_name"],
        pbmm2_threads=config["pbmm2_threads"]
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 pbmm2 align --preset ISOSEQ -O6,24 -B4 --sort /data/{input.isoseq3_reads} /data/{input.bedtools_output} /data/{output.pbmm2_output} -j {params.pbmm2_threads}"

# Collapse redundant transcrips from the pbmm2 outputs
rule isoseq3:
    input:
        pbmm2_output=rules.pbmm2.output
    output:
        isoseq3_output=f"isoseq_outputs/{species_name}.collapsed.gff"
    params:
        volume_name=config["volume_name"],
        isoseq3_min_coverage=config["isoseq3_min_coverage"],
        isoseq3_min_identity=config["isoseq3_min_identity"]
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 isoseq3 collapse --do-not-collapse-extra-5exons --min-aln-coverage {params.isoseq3_min_coverage} --min-aln-identity {params.isoseq3_min_identity} /data/{input.pbmm2_output} /data/{output}"

# Create maker control files which are the inputs for maker
rule maker_control_files:
    input:
        isoseq3_output=f"isoseq_outputs/{species_name}.collapsed.gff"
    output:
        chrs_list= f"maker_outputs/chrs.list"
    params:
        fasta_file=f"inputs/{config['fasta_file']}",
        isoseq3_output=f"isoseq_outputs/{species_name}.collapsed.gff",
        protein_file=f"inputs/{config['protein_file']}",
        repeatmodeler_output=f"repeat_modeler_outputs/{species_name}-families.fa",
        volume_name=config["volume_name"]
    shell:
        "bash {params.volume_name}/scripts/maker_control_files {params.fasta_file} {params.isoseq3_output} {params.protein_file} {params.repeatmodeler_output} {params.volume_name} > {output.chrs_list}"

# Create genome annotations 
rule maker:
    input:
        done = rules.maker_control_files.output
    output:
        done = touch("maker_outputs/{chr}/maker_done")
    shell:
        "docker run -v {config[volume_name]}:/data danifilho/danifilho:latest /bin/bash -c 'cd /data/maker_outputs/{wildcards.chr} && maker /data/maker_outputs/{wildcards.chr}/maker_opts.ctl /data/maker_outputs/{wildcards.chr}/maker_bopts.ctl /data/maker_outputs/{wildcards.chr}/maker_exe.ctl'"


# Creates a global GFF file by concatenating all the individual GFF files generated by MAKER for each chromosome. Filters for only gene, CDS, mRNA, exon, UTR, and tRNA features.
rule create_global_gff:
    input:
        maker_done = expand("maker_outputs/{chr}/maker_done", chr=chr)
    output:
        global_gff =f"filtering_outputs/{config['species_abbreviation']}.v01.gff"
    shell:
        """
        printf "##gff-version 3\n" > {output}
        tail -n +2 maker_outputs/*/*/*_datastore/*/*/*/*.gff | awk -F"\t" 'NF==9 && ($3=="gene" || $3 =="CDS" || $3 =="mRNA" || $3 =="exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" || $3 =="tRNA" )' >> {output}
        """

# Extracts protein sequences from the global GFF file using gffread. Also indexes the protein FASTA file with samtools.
rule extract_proteins_and_index:
    input:
        global_gff = rules.create_global_gff.output
    output:
        prots_fai = f"filtering_outputs/{config['species_abbreviation']}.prots.fai"
    params:
        fasta_file=f"inputs/{config['fasta_file']}",
        volume_name = config['volume_name'],
        species_abbreviation = config['species_abbreviation']
    shell:
        """
        docker run -v {params.volume_name}:/data biodepot/gffread gffread -S -y "/data/filtering_outputs/{params.species_abbreviation}.prots" -g "/data/{params.fasta_file}" "/data/{input.global_gff}"
        docker run -v {params.volume_name}:/data dbest/samtools:v1.19.2 samtools faidx "/data/filtering_outputs/{params.species_abbreviation}.prots"
        """

# Finds the smallest protein length value in the protein FASTA index. Uses this value to generate a list of protein IDs below this threshold to exclude.
rule find_and_process_smallest_value:
    input:
        prots_fai = rules.extract_proteins_and_index.output
    output:
        exclude_list= f"filtering_outputs/exclude.list"
    shell:
        """
        sort -k2 -n {input.prots_fai} > sorted_prots.fai
        MENOR_VALOR=$(sort -r -k2 -n {input.prots_fai} | cut -f -2 | grep "est2genome" | tail -n 1 | cut -f2)
        sort -k2 -n {input.prots_fai} | cut -f -2 | grep protein2genome | awk -v minor="$MENOR_VALOR" '{{if ($2 < minor) print $1}}' | sed 's/-mRNA-[0-9]*//' | sort | uniq > {output.exclude_list}
        rm sorted_prots.fai
        """

# Filters the original GFF to remove proteins below length threshold. Validates and makes IDs unique. Renames IDs with species abbreviation prefix.
rule process_script:
    input:
        prots_fai = f"filtering_outputs/{config['species_abbreviation']}.prots.fai",
        exclude_list= rules.find_and_process_smallest_value.output
    output:
        exclude_gff = f"filtering_outputs/{config['species_abbreviation']}.exclude.gff",
        uniq_gff = f"filtering_outputs/{config['species_abbreviation']}.uniq.gff",
        final_gff = f"filtering_outputs/{config['species_abbreviation']}.v2.gff"
    params:
        volume_name = config["volume_name"],
        global_gff = f"filtering_outputs/{config['species_abbreviation']}.v01.gff",
        species_abbreviation = config['species_abbreviation']
    shell:
        """
        printf "#gff-version 3\n" > {output.exclude_gff}
        docker run -v {params.volume_name}:/data python sh -c "pip install gff3 && python /data/scripts/gff_filter.py -e /data/{input.exclude_list} -g /data/{params.global_gff} >> /data/{output.exclude_gff}"
        docker run -v {params.volume_name}:/data dantestpy1 sh -c "python /data/scripts/validate_gff.py --gff /data/{output.exclude_gff} > {output.uniq_gff}"
        docker run -v {params.volume_name}:/data dantestpy2 sh -c "python /data/scripts/renameGff.py -g /data/{output.uniq_gff} -t /data/{params.species_abbreviation} > {output.final_gff}"
        """
