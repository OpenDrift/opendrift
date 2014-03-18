from readers import PointReader

class Reader(PointReader):

    # Parameters (CF standard names) which
    # can be provided by this model/reader
    parameters = ['eastward_wind',
                  'northward_wind',
                  'air_temperature']
        
    def __init__(self):
        # e.g. open file, check that everything is ok
        pass
        
    def read_parameters(self, time=None, lon=None, lat=None,
                        depth=None, parameters):
        # Fake eastward wind with velocity 10 m/s
        if 'eastward_wind' in parameters:
            self.easthward_wind = 10
        if 'northward_wind' in parameters:
            self.northward_wind = 0
