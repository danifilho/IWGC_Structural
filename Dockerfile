# Usando a imagem base do alpine/git para clonar o repositório
FROM alpine/git AS base-alpine

# Definindo o diretório de trabalho
WORKDIR /myrepo

# Instalando o 'make' e clonando o repositório Git
RUN apk add --update make && \
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && \
    make augustus
    
# Usando a imagem base do miniconda 
FROM continuumio/miniconda3

# Criando um diretório de trabalho
WORKDIR /app

# Copiando o repositório clonado da imagem base-alpine
COPY --from=base-alpine /myrepo /app/myrepo

# Copiando o arquivo environment.yml para o contêiner
COPY environment.yml .

# Usando o conda para instalar o Maker
RUN conda env create -f environment.yml && conda clean -a

# Ativando o ambiente conda
RUN conda init bash

# Definindo o comando para iniciar o Maker
CMD ["maker"]
