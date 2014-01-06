# Prototype opendrift module, containing two generic classes:
#   - LagrangianElement
#   - OpenDriftSimulation
#
# Specific modules (e.g. oildrift or leeway) will inherit from these classes
#
# Knut-Frode Dagestad, Jan 2014


class LagrangianElement():
    ''' Class which stores data for each LagrangianElement '''

    def __init__(self, lon=None, lat=None):
        print 'Initialising a generic LagrangianElement'
        self.lon = lon
        self.lat = lat


class OpenDriftSimulation():
    ''' Class which stores data (e.g. lagrangian elements) and methods for a drift simulation '''
    
    def __init__(self):
        print 'Initialising a generic drift simulation'
        self.activeLagrangianElements = []

    def add_element(self, lon, lat):
        self.activeLagrangianElements.append(LagrangianElement(lon, lat))

    def print_element_positions(self):
        print 'Reporting element positions:'
        for i, activeElement in enumerate(self.activeLagrangianElements):
            print '\tElement number %i at position (%f, %f)' % (i, activeElement.lon, activeElement.lat)

    def run_drift_simulation(self):
        print '\nRunning a drift simulation...'
        self.print_element_positions()
        #
        # Loop:
        # - seed new elements
        # - update particle properties (e.g. oil weathering or leeway jibing) 
        # - advect elements
        # - deactivate stranded elements
        # - check end conditions
        #   - no more active elements
        #   - end of driver data period
