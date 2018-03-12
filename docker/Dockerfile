FROM continuumio/miniconda
MAINTAINER Simon Weppe <simon@metocean.co.nz>

RUN apt install vim libgl1-mesa-glx

# Install Opendrift
# This is using my forked version with the edits on the readers/examples for now
RUN echo "-----------------Install OPENDRIFT-----------------------------------------------" &&\
    conda install --yes hdf4 numpy scipy matplotlib basemap netcdf4 configobj pillow gdal pyproj &&\
    git clone https://github.com/simon-msl/opendrift.git /source/opendrift  &&\
    cd /source/opendrift &&\
    python setup.py develop --user &&\
    export PATH=$PATH:$OPENDRIFT_FOLDER/opendrift/scripts/



# Install NOAA's OilLibrary
# requires setuptools and some other packages
# many packages available on anaconda's noaa-orr-erd channel 
RUN echo "-----------------Install NOAA's OilLibrary----------------------------------------" &&\
    git clone https://github.com/NOAA-ORR-ERD/OilLibrary.git /source/OilLibrary &&\
    wget https://bootstrap.pypa.io/ez_setup.py -O - | python &&\ 
    conda install --yes sqlalchemy transaction zope.sqlalchemy pytest &&\
    conda install --yes -c anaconda repoze.lru &&\
    conda install --yes -c noaa-orr-erd awesome-slugify  unit_conversion &&\
    cd /source/OilLibrary &&\
    python setup.py install


# Install motu-client to be able to download MERCATOR/CMEMS data
# https://github.com/clstoulouse/motu-client-python#Configuration
RUN echo "-----------------Install Motu-client - For Mercator Data downloads-----------------" &&\
    cd /source &&\
    wget "https://github.com/clstoulouse/motu-client-python/releases/download/motu-client-python-1.5.00/motu-client-python-1.5.00-20180223190259664-bin.tar.gz" &&\
    tar xvzf motu-client-python-1.5.00-20180223190259664-bin.tar.gz &&\
    rm *.gz