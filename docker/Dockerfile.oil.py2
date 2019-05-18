FROM continuumio/miniconda

MAINTAINER vsochat@stanford.edu

# https://hub.docker.com/r/vanessa/opendrift-oil/  # container 
# https://github.com/opendrift/opendrift/wiki      # documentation

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
RUN /opt/conda/bin/conda env create -f conda_python2_oil.yml
RUN /bin/bash -c '''. activate opendrift_p2_oil && \
    python setup.py develop && \
    echo "source activate opendrift_p2_oil" > ~/.bashrc'''

# NOAA's OilLibrary for Oil Simulations
WORKDIR /
RUN git clone https://github.com/NOAA-ORR-ERD/OilLibrary.git && \
    /opt/conda/bin/conda install --yes sqlalchemy transaction zope.sqlalchemy pytest && \
    /opt/conda/bin/conda install --yes -c anaconda repoze.lru && \
    /opt/conda/bin/conda install --yes -c noaa-orr-erd awesome-slugify unit_conversion && \
    cd OilLibrary && \
    /bin/bash -c '''. activate opendrift_p2_oil && \
    pip install gitpython && \
    python setup.py install'''

WORKDIR /code

# Test
# RUN /bin/bash -c ". activate opendrift_p2_oil && cd /code && ./testall"
