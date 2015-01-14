from readers import readers

class OpenDriftSimulation(object):

    def __init__(self):


        # Add Readers object as single point interface 
        # towards environmental data
        self.readers = readers.Readers()

        print 'OpenDriftSimulation initialised'
