#!/usr/bin/env python
"""
Fetch positions of drifting buoys and simulate trajectory
=========================================================
"""


#%%
# Requires sfy-processing from sfy
#
# pip install "git+https://github.com/gauteh/sfy#egg=sfy-processing&subdirectory=sfy-processing"

import sfy
import cmocean
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

hub = sfy.hub.Hub.from_env()

