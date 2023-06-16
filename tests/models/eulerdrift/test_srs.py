import pytest
import utm

from opendrift.models import eulerdrift

def test_norway():
  p = eulerdrift.srs.find_utm_proj(9, 60)
  assert 'zone=32' in str(p)

def test_oor():
  with pytest.raises(utm.error.OutOfRangeError):
    _ = eulerdrift.srs.find_utm_proj(0, 90)
