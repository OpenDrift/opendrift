# See https://opendrift.github.io for usage

FROM condaforge/mambaforge

ENV DEBIAN_FRONTEND noninteractive

RUN mkdir /code
WORKDIR /code

# Install opendrift environment into base conda environment
COPY environment.yml .
RUN mamba env update -n base -f environment.yml

# Cache cartopy maps
RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Install opendrift
ADD . /code
RUN pip install -e .

# Test installation
RUN /bin/bash -c "echo -e \"import opendrift\" | python"

