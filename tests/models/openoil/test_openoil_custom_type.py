from datetime import datetime
import json
from importlib.resources import files

from opendrift.models.openoil import OpenOil

oiljson = files('opendrift.models.openoil.adios.extra_oils').joinpath('AD04011.json')

def test_set_oil_type_json():
    o = OpenOil(loglevel=50)

    with open(oiljson) as f:
        d = json.loads(f.read())

    o.set_oiltype_by_json(d)
    assert o.oiltype.name == 'GENERIC FUEL OIL No. 6'
    print(o.oiltype)

    o = OpenOil(loglevel=50)
    o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0,
                                          'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0})
    o.seed_elements(lon=0, lat=0, time=datetime.now(), oil_type=d)
    assert o.oiltype.name == 'GENERIC FUEL OIL No. 6'

def test_set_oil_type_file():
    o = OpenOil(loglevel=50)

    o.set_oiltype_from_file(oiljson)

    assert o.oiltype.name == 'GENERIC FUEL OIL No. 6'
    print(o.oiltype)

    o = OpenOil(loglevel=50)
    o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0,
                                          'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0})
    o.seed_elements(lon=0, lat=0, time=datetime.now(), oil_type=oiljson)
    assert o.oiltype.name == 'GENERIC FUEL OIL No. 6'
