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
ADD . .
RUN micromamba run -n base pip install -e .


# Download cartopy GSHHS shoreline - temporary solution until Cartopy is updated
USER root
RUN /bin/bash -c "mkdir -p /var/lib/apt/lists/partial && \
    DEBIAN_FRONTEND=noninteractive apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y unzip wget && \
    mkdir -p /root/.local/share/cartopy/shapefiles/ && \
    cd /root/.local/share/cartopy/shapefiles/ && \
    rm -f gshhg-shp-2.3.7.zip && \
    rm -rf gshhs && \
    wget https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip && \
    unzip -o gshhg-shp-2.3.7.zip && \
    mv GSHHS_shp gshhs"
USER $MAMBA_USER

## Cache cartopy maps
#RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Test installation
RUN micromamba run -n base python -c "import opendrift"
