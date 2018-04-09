****WORK IN PROGRESS ***

INSTALL NOTES 
------------------

# to build the docker image :

cd /docker
docker build -t metocean/opendrift:simon .

run (interactively) : 

docker run -it metocean/opendrift:simon
docker run -it -v /home/simon:/home/simon metocean/opendrift:simon # example with a mounted folder

------------------
Following instruction from :
https://github.com/opendrift/opendrift/wiki

git clone https://github.com/OpenDrift/opendrift.git
cd opendrift; 
python setup.py develop --user
export PATH=$PATH:$OPENDRIFT_FOLDER/opendrift/scripts/
conda install --yes hdf4 numpy scipy matplotlib basemap netcdf4 configobj pillow gdal pyproj