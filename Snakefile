# Snakemake file: Snakefile
configfile: "config.yaml"

from Bio import SeqIO
chr = []
for record in SeqIO.parse(config["fasta_file"], "fasta"):
    chr.append(record.id)
    
rule all:
    input:
        expand(f"{config['species_abbreviation']}.v2.gff")

rule repeatmodeler:
    input:
        fasta_file=config["fasta_file"]
    output:
        repeatmodeler_output=f"{config['species_name']}-families.fa"
    params:
        volume_name=config["volume_name"],
        num_jobs_repeat=config["num_jobs_repeat"],
        database_name=f"{config['species_name']}"
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest BuildDatabase -name /data/{params.database_name} /data/{input.fasta_file} && "
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatModeler -database /data/{params.database_name} -threads {params.num_jobs_repeat} -LTRStruct"

rule repeatmasker:
    input:
        fasta_file=config["fasta_file"],
        repeatmodeler_output=rules.repeatmodeler.output
    output:
        repeatmasker_output=f"{config['species_name']}.fasta.out.gff"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"]
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatMasker -gff -a -pa 16 -u /data/{input.repeatmodeler_output} /data/{input.fasta_file}"

rule bedtools:
    input:
        fasta_file=config["fasta_file"],
        repeatmasker_output=rules.repeatmasker.output
    output:
        bedtools_output=f"{config['species_name']}.softmask.fasta"
    params:
        volume_name=config["volume_name"]
    shell:
        "docker run -v \"{params.volume_name}\":/data pegi3s/bedtools:latest bedtools maskfasta -fi /data/{input.fasta_file} -bed /data/{input.repeatmasker_output} -soft -fo /data/{output.bedtools_output}"

rule pbmm2:
    input:
        bedtools_output=rules.bedtools.output
    output:
        pbmm2_output=f"{config['species_name']}-aligned.bam"
    params:
        volume_name=config["volume_name"],
        isoseq3_reads=config["isoseq3_reads"],
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 pbmm2 align --preset ISOSEQ -O6,24 -B4 --sort /data/{params.isoseq3_reads} /data/{input.bedtools_output} /data/{output.pbmm2_output} -j 8"

rule isoseq3:
    input:
        pbmm2_output=rules.pbmm2.output
    output:
        isoseq3_output=f"{config['species_name']}.collapsed.gff"
    params:
        volume_name=config["volume_name"],
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 isoseq3 collapse --do-not-collapse-extra-5exons --min-aln-coverage 0.9 --min-aln-identity 0.95 /data/{input.pbmm2_output} /data/{output}"

rule maker_control_files:
    input:
        done=rules.isoseq3.output
    output:
        "chrs.list"
    params:
        fasta_file=config["fasta_file"],
        isoseq3_output=f"{config['species_name']}.collapsed.gff",
        protein_file=config["protein_file"],
        repeatmodeler_output=f"{config['species_name']}-families.fa",
        volume_name=config["volume_name"]
    shell:
        "bash maker_control_files {params.fasta_file} {params.isoseq3_output} {params.protein_file} {params.repeatmodeler_output} {params.volume_name}> chrs.list"

rule maker:
    input:
        done = rules.maker_control_files.output
    output:
        done = touch("{volume_name}/{chr}/maker_done")
    shell:
        "docker run -v {config[volume_name]}:/data danifilho/danifilho:latest /bin/bash -c 'cd /data/{wildcards.chr} && maker /data/{wildcards.chr}/maker_opts.ctl /data/{wildcards.chr}/maker_bopts.ctl /data/{wildcards.chr}/maker_exe.ctl'"

rule create_global_gff:
    input:
        expand("{volume_name}/{chr}/maker_done", 
        chr=chr, 
        volume_name=config["volume_name"])
    output:
        global_gff =f"{config['species_abbreviation']}.v01.gff"
    shell:
        """
        printf "##gff-version 3\n" > {output}
        tail -n +2 */*maker.output/*_datastore/*/*/*/*gff | awk -F"\t" 'NF==9 && ($3=="gene" || $3 =="CDS" || $3 =="mRNA" || $3 =="exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" || $3 =="tRNA" )' >> {output}
        """

rule extract_proteins_and_index:
    input:
        global_gff = rules.create_global_gff.output
    output:
        prots_fai = f"{config['species_abbreviation']}.prots.fai"
    params:
        fasta_file=config["fasta_file"],
        volume_name = config['volume_name'],
        species_abbreviation = config['species_abbreviation']
    shell:
        """
        docker run -v {params.volume_name}:/data biodepot/gffread gffread -S -y "/data/{params.species_abbreviation}.prots" -g "/data/{params.fasta_file}" "/data/{input.global_gff}"
        docker run -v {params.volume_name}:/data dbest/samtools:v1.19.2 samtools faidx "/data/{params.species_abbreviation}.prots"
        """

rule find_and_process_smallest_value:
    input:
        prots_fai = rules.extract_proteins_and_index.output
    output:
        "exclude.list"
    shell:
        """
        sort -k2 -n {input.prots_fai} > sorted_prots.fai
        MENOR_VALOR=$(sort -r -k2 -n {input.prots_fai} | cut -f -2 | grep "est2genome" | tail -n 1 | cut -f2)
        sort -k2 -n {input.prots_fai} | cut -f -2 | grep protein2genome | awk -v minor="$MENOR_VALOR" '{{if ($2 < minor) print $1}}' | sed 's/-mRNA-[0-9]*//' | sort | uniq > exclude.list
        rm sorted_prots.fai
        """

rule process_script:
    input:
        prots_fai = f"{config['species_abbreviation']}.prots.fai",
        exclude_list="exclude.list"
    output:
        exclude_gff = f"{config['species_abbreviation']}.exclude.gff",
        uniq_gff = f"{config['species_abbreviation']}.uniq.gff",
        final_gff = f"{config['species_abbreviation']}.v2.gff"
    params:
        volume_name = config["volume_name"],
        global_gff = f"{config['species_abbreviation']}.v01.gff",
        species_abbreviation = config['species_abbreviation']
    shell:
        """
        printf "#gff-version 3\n" > {output.exclude_gff}
        docker run -v {params.volume_name}:/data python sh -c "pip install gff3 && python /data/gff_filter.py -e /data/{input.exclude_list} -g /data/{params.global_gff} >> /data/{output.exclude_gff}"
        docker run -v {params.volume_name}:/data dantestpy1 sh -c "python /data/validate_gff.py --gff /data/{output.exclude_gff} > {output.uniq_gff}"
        docker run -v {params.volume_name}:/data dantestpy2 sh -c "python /data/renameGff.py -g /data/{output.uniq_gff} -t /data/{params.species_abbreviation} > {output.final_gff}"
        """
