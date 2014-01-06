#!/usr/bin/env python
# Test script to illustrate usage
# Knut-Frode Dagestad, Jan 2014

from opendrift import *

# Initialising a new drift simulation object
os = OpenDriftSimulation()

# Add two LagrangianElements, outside Bergen
os.add_element(5, 60)
os.add_element(5.1, 60.1)

# Finally, run simulation (presently only reporting locations of the elements)
os.run_drift_simulation()
