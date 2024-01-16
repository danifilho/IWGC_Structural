# Use a imagem base do Ubuntu
FROM ubuntu:latest

# Atualize o sistema e instale as dependências necessárias
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git build-essential wget cpanminus zlib1g-dev libexpat1-dev

# Clone o repositório do Maker
RUN git clone https://github.com/Yandell-Lab/maker.git

# Mude para o diretório /src do Maker
WORKDIR /maker/src/

# Instale o Maker
RUN cpanm Thread::Queue && \
    perl Build.PL && \
    ./Build installdeps && \
    ./Build install

ENV PATH /maker/bin:$PATH
