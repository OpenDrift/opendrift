from readers import PointReader

class Reader(PointReader):

    # Parameters (CF standard names) which
    # can be provided by this model/reader
    parameters = ['eastward_sea_water_velocity',
                  'northward_sea_water_velocity',
                  'sea_water_salinity',
                  'sea_water_temperature']

    def __init__(self):
        # e.g. open file, check that everything is ok
        pass

    def read_parameters(self, parameters, time=None, lon=None, lat=None,
                        depth=None ):
        # Fake northward current with velocity 1 m/s
        if 'eastward_sea_water_velocity' in parameters:
            self.easthward_sea_water_velocity = 0
        if 'northward_sea_water_velocity' in parameters:
            self.northward_sea_water_velocity = 1
