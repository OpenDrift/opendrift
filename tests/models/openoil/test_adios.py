from opendrift.models.openoil import adios

def test_get_all():
    oils = adios.oils(300)
    assert len(oils) == 300

    # make sure they are all unique
    ids = [o._id for o in oils]
    sids = set(ids)
    assert len(ids) == len(sids)
