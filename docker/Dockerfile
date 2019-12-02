FROM continuumio/miniconda3

MAINTAINER vsochat@stanford.edu

## The container is provided on Docker Hub
# https://hub.docker.com/r/opendrift/opendrift

## Documentation for opendrift is here
# https://github.com/opendrift/opendrift/wiki

ENV DEBIAN_FRONTEND noninteractive
ENV PATH /code/opendrift/opendrift/scripts:$PATH

# Dependencies for opendrift
RUN apt-get update && apt-get install -y build-essential \
                                         apt-utils \
                                         unzip \
                                         vim \
                                         git \
                                         gfortran \
                                         libgeos-dev \
                                         gdal-bin && \
                                         ldconfig

RUN pip install --upgrade pip
RUN mkdir /code
ADD . /code
WORKDIR /code
RUN /opt/conda/bin/conda env create -f environment.yml
RUN /bin/bash -c '''. activate opendrift && \
    pip install -e . && \
    echo "source activate opendrift" > ~/.bashrc'''

WORKDIR /code

# Test
# RUN /bin/bash -c ". activate opendrift && cd /code && ./testall"
