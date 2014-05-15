from readers import PointReader
from netCDF4 import Dataset

class Reader(PointReader):

    # Parameters (CF standard names) which
    # can be provided by this model/reader
    parameters = ['x_sea_water_velocity',
                  'y_sea_water_velocity',
                  'sea_water_salinity',
                  'sea_water_temperature']

    def __init__(self):
        # hardcoded filename
        f = 'http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an.2014032500.nc'
        # Open file, check that everything is ok
        self.Dataset = Dataset(f)
        # Read and store arrays of longitude and latitude
        self.lon = self.Dataset.variables['lon'][:]
        self.lat = self.Dataset.variables['lat'][:]
        print self.lat.shape


    def read_parameters(self, parameters, time=None, lon=None, lat=None,
                        depth=None ):
        # Fake northward current with velocity 1 m/s
        if 'eastward_sea_water_velocity' in parameters:
            self.easthward_sea_water_velocity = 0
        if 'northward_sea_water_velocity' in parameters:
            self.northward_sea_water_velocity = 1