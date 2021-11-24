#! /usr/bin/env python
"""
Example plot (for front page)
==================================
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=30)
o.add_readers_from_list(
    ['https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])
o.disable_vertical_motion()
o.seed_elements(lon=4.85, lat=60, time=datetime.now(), number=10000, radius=1000)

o.run(duration=timedelta(hours=24))
o.animation(filename='animation.mp4')

#%%
# .. image:: /gallery/animations/example_plot_0.gif

