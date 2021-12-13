# See https://opendrift.github.io for usage

FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND noninteractive
ENV PATH /code/opendrift/opendrift/scripts:$PATH

RUN mkdir /code
WORKDIR /code

RUN conda config --add channels noaa-orr-erd
RUN conda config --add channels conda-forge
RUN conda config --add channels opendrift

# Install opendrift environment into base conda environment
COPY environment.yml .
RUN /opt/conda/bin/conda env update -n base -f environment.yml

# Install roaring-landmask
RUN pip install roaring-landmask

# Cache cartopy maps
RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Install opendrift
ADD . /code
RUN pip install -e .

# Test installation
RUN /bin/bash -c "echo -e \"import opendrift\" | python"

