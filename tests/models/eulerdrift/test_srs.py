import pytest
import advent
import utm

def test_norway():
  p = advent.srs.find_utm_proj(9, 60)
  assert 'zone=32' in str(p)

def test_oor():
  with pytest.raises(utm.error.OutOfRangeError):
    _ = advent.srs.find_utm_proj(0, 90)
