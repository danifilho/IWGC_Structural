# Snakemake file: Snakefile
configfile: "config.yaml"

rule all:
    input:
        expand("{volume_name}/{species_name}-families.fa", volume_name=config["volume_name"], species_name=config["species_name"])

rule repeatmodeler:
    input:
        fasta_file=config["fasta_file"]
    output:
        repeatmodeler_output="{volume_name}/{species_name}-families.fa"
    params:
        species_name=config["species_name"],
        volume_name=config["volume_name"],
        database_name=config["database_name"],
        num_jobs_repeat=config["num_jobs_repeat"]
    shell:
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest BuildDatabase -name /data/{params.database_name} /data/{input.fasta_file} && "
        "docker run -v \"{params.volume_name}\":/data dfam/tetools:latest RepeatModeler -database /data/{params.database_name} -threads {params.num_jobs_repeat} -LTRStruct"
