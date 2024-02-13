# Snakemake file: Snakefile
configfile: "config.yaml"

rule all:
    input:
       expand("{volume_name}/{species_name}.collapsed.gff", volume_name=config["volume_name"], species_name=config["species_name"])
rule repeatmodeler:
    input:
        fasta_file=config["fasta_file"]
    output:
        repeatmodeler_output=config["repeatmodeler_output"]
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"],
        database_name=config["database_name"],
        num_jobs_repeat=config["num_jobs_repeat"]
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest BuildDatabase -name /data/{params.database_name} /data/{input.fasta_file} && "
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatModeler -database /data/{params.database_name} -threads {params.num_jobs_repeat} -LTRStruct"

rule repeatmasker:
    input:
        fasta_file=config["fasta_file"],
        repeatmodeler_output=config["repeatmodeler_output"]

    output:
        repeatmasker_output="{volume_name}/{species_name}.fasta.out.gff"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"]
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatMasker -gff -a -pa 16 -u /data/{input.repeatmodeler_output} /data/{input.fasta_file}"

rule bedtools:
    input:
        fasta_file="file1.Col-PEK1.5_Chr1-5_20220523.fasta",
        repeatmasker_output="file1.Col-PEK1.5_Chr1-5_20220523.fasta.out.gff"
    output:
        bedtools_output="{volume_name}/{species_name}.softmask.fasta"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"],
        bedtools_output=config["bedtools_output"]
    shell:
        "docker run -v \"{params.volume_name}\":/data pegi3s/bedtools:latest bedtools maskfasta -fi /data/{input.fasta_file} -bed /data/{input.repeatmasker_output} -soft -fo /data/{params.bedtools_output}"

rule pbmm2:
    input:
        isoseq3_reads=config["isoseq3_reads"],
        bedtools_output=config["bedtools_output"]

    output:
        pbmm2_output="{volume_name}/{species_name}-aligned.bam"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"],
        pbmm2_output=config["pbmm2_output"]
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 pbmm2 align --preset ISOSEQ -O6,24 -B4 --sort /data/{input.isoseq3_reads} /data/{input.bedtools_output} /data/{params.pbmm2_output} -j 8"

rule isoseq3:
    input:
        pbmm2_output=config["pbmm2_output"]

    output:
        isoseq3_output="{volume_name}/{species_name}.collapsed.gff"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"],
        isoseq3_output=config["isoseq3_output"]
    shell:
        "docker run -v \"{params.volume_name}\":/data greensii/isoseq3 isoseq3 collapse --do-not-collapse-extra-5exons --min-aln-coverage 0.9 --min-aln-identity 0.95  /data/{input.pbmm2_output} /data/{params.isoseq3_output}"

#---------------- filtering working

rule process_script:
    input:
        prots_fai = f"{species_abbreviation}.prots.fai",
        exclude_list = "exclude.list"
    output:
        exclude_gff = f"{species_abbreviation}.exclude.gff",
        uniq_gff = f"{species_abbreviation}.uniq.gff",
        final_gff = f"{species_abbreviation}.v2.gff"
    params:
        volume_name = volume_name
    shell:
        """
        printf "#gff-version 3\n" > {output.exclude_gff}
        docker run -v {params.volume_name}:/data python sh -c "pip install gff3 && python /data/gff_filter.py -e /data/{input.exclude_list} -g /data/{global_gff} >> /data/{output.exclude_gff}"
        docker run -v {params.volume_name}:/data dantestpy sh -c "python /data/validate_gff.py --gff /data/{output.exclude_gff} > {output.uniq_gff}"
        docker run -v {params.volume_name}:/data dantestpy sh -c "python /data/renameGff.py -g /data/{output.uniq_gff} -t /data/{species_abbreviation} > {output.final_gff}"
        """
