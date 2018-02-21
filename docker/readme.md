INSTALL NOTES 
------------------

As root build using 
/home/simon/docker/ercore
docker build -t metocean/opendrift_simon .

run (interactively)

docker run -it metocean/opendrift_simon

docker run -it -v /home/simon:/home/simon metocean/opendrift_simon # example with a mounted folder


------------------
https://github.com/opendrift/opendrift/wiki

git clone https://github.com/OpenDrift/opendrift.git
cd opendrift; 
python setup.py develop --user
export PATH=$PATH:$OPENDRIFT_FOLDER/opendrift/scripts/
conda install --yes hdf4 numpy scipy matplotlib basemap netcdf4 configobj pillow gdal pyproj

conda is the package manager. Anaconda is a set of about a hundred packages including conda, numpy, scipy, ipython notebook, and so on. 