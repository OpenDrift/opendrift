from opendrift.models import eulerdrift

def test_contains():
  grid = eulerdrift.grid.RegularGrid.new(10, 65, 10, (1000, 1000))
  x, y = grid.srs(10.1, 65.01)

  assert grid.contains(x, y)
