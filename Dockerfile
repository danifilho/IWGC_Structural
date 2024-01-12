# Usando a imagem base do Ubuntu 
FROM ubuntu:latest AS base-ubuntu

WORKDIR /myrepo

# Atualizando a lista de pacotes e instalando o 'git' e o 'make'
RUN apt-get update && apt-get install -y git make build-essential wget git autoconf libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev libbamtools-dev zlib1g-dev && \
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && \
    make augustus

# Usando a imagem base do miniconda 
FROM continuumio/miniconda3

# Criando um diretório de trabalho
WORKDIR /app

# Copiando o repositório clonado da imagem base-ubuntu
COPY --from=base-ubuntu /myrepo/Augustus/bin/augustus /usr/local/bin/augustus

# Copiando o arquivo environment.yml para o contêiner
COPY environment.yml .

# Usando o conda para instalar o Maker
RUN conda env create -f environment.yml && conda clean -a

# Ativando o ambiente conda
RUN conda init bash

# Definindo o comando para iniciar o Maker
CMD ["maker"]
