# See https://opendrift.github.io for usage

# Use a minimal base image
FROM mambaorg/micromamba:2.1.1

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_DOCKERFILE_ACTIVATE=1

RUN mkdir code
WORKDIR code

# Copy environment file
COPY environment.yml .

# Install opendrift environment into base micromamba environment
RUN micromamba install -n base -f environment.yml

# Cache cartopy maps
RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Install opendrift
ADD . .
RUN pip install -e .

# Test installation
RUN /bin/bash -c "echo -e \"import opendrift\" | python"
