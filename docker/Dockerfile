FROM continuumio/miniconda
MAINTAINER Simon Weppe <simon@metocean.co.nz>


RUN echo "-----------------Install OPENDRIFT-----------------" &&\
    conda install --yes hdf4 numpy scipy matplotlib basemap netcdf4 configobj pillow gdal pyproj &&\
    git clone https://github.com/simon-msl/opendrift.git /opendrift  &&\
    cd /opendrift &&\
    python setup.py develop --user &&\
    export PATH=$PATH:$OPENDRIFT_FOLDER/opendrift/scripts/

# This is using my forked version with the edits on the readers/examples etc..
# 
#
# some issues witg GDAL - says it is not recognized ?
#
# use this to use own branch : git clone https://github.com/simon-msl/opendrift.git /opendrift
#
#
# 
# Would be useful to install the motu_client.py
# available from https://github.com/clstoulouse/motu-client-python
# To download data directly from CMEMS
#
#
# Install NOAA's Oil libray as well


