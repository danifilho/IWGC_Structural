# Base Image
FROM continuumio/anaconda3

# Software Version
LABEL maker="2.31.10"

# Description
LABEL description="MAKER is a portable and easily configurable genome annotation pipeline."

# Update conda
RUN conda update -n base -c defaults conda -y

# Create a new conda environment with Python 2.7
RUN conda create -n maker_env python=2.7 -y

# Activate the new environment and install maker
RUN echo "source activate maker_env" > ~/.bashrc
ENV PATH /opt/conda/envs/maker_env/bin:$PATH

# Install h5py
RUN /bin/bash -c "source activate maker_env && conda config --add channels conda-forge"
RUN /bin/bash -c "source activate maker_env && conda install -c anaconda h5py -y"
RUN /bin/bash -c "source activate maker_env && conda install -c anaconda gcc_linux-64 -y"
RUN /bin/bash -c "source activate maker_env && conda install -c anaconda gfortran_linux-64 -y"
RUN /bin/bash -c "source activate maker_env && conda install -c anaconda gxx_linux-64 -y"

# Install maker
RUN /bin/bash -c "source activate maker_env && conda install -c bioconda maker -y"

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD . /app
