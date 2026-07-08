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

# Download cartopy GSHHS shoreline - temporary solution until Cartopy is updated
SHELL ["/bin/bash", "-c"]
USER root
RUN DEBIAN_FRONTEND=noninteractive apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y unzip wget && \
    source activate opendrift && \
    mkdir -p /root/.local/share/cartopy/shapefiles/ && \
    cd /root/.local/share/cartopy/shapefiles/ && \
    rm -f gshhg-shp-2.3.7.zip && \
    rm -rf gshhs && \
    wget https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip && \
    unzip -o gshhg-shp-2.3.7.zip && \
    mv GSHHS_shp gshhs

## Cache cartopy maps
#RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Install opendrift
ADD . .
RUN pip install -e .

# Test installation
RUN /bin/bash -c "echo -e \"import opendrift\" | python"
