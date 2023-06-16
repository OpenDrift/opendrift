import numpy as np
from opendrift.models import eulerdrift

def test_init():
  s = eulerdrift.Simulation.new()
  s.grid.grid[0:10,0:10] = 1.

def test_source():
  s = eulerdrift.Simulation.new()
  s.source(s.grid.lon0, s.grid.lat0, np.ones((10, 10)))
  assert s.grid.grid[0,0] == 1.

def test_source_gaussian():
  s = eulerdrift.ExplSimulation.new()
  s.readers.append(eulerdrift.ConstantReader.new_xy())

  loc, lac = s.grid.center()

  s.source_gaussian_blob(loc, lac, 1., 70, 20.)

  s.grid.plot()
  s.integrate(dt = .01, max_steps = 5000)

  s.grid.plot()
