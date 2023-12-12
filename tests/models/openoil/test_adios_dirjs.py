import pytest
from opendrift.models.openoil import adios


def test_get_oils(benchmark):
    oils = benchmark(adios.dirjs.oils)
    assert len(oils) == 50


def test_get_all_oils(benchmark):
    oils = benchmark(adios.dirjs.oils, None)
    assert len(oils) >= 1269


def test_get_all_oil_names(benchmark):
    oils = benchmark(adios.dirjs.get_oil_names)
    assert len(oils) >= 1269


def test_get_all_oil_names_location(benchmark):
    oils = benchmark(adios.dirjs.get_oil_names, location='NORWAY')
    assert len(oils) >= 150 and len(oils) <= 201


def test_bunker_c_1987():
    o = adios.dirjs.oils(query='Bunker C [1987]')
    print([oo.name for oo in o])
    assert len(o) == 1
    assert o[0].valid()


def test_get_oil_by_name(benchmark):
    oil = benchmark(adios.dirjs.find_full_oil_from_name, 'AASGARD A 2003')
    assert oil.name == 'AASGARD A 2003'
    assert oil.valid()


def test_get_oil_by_id(benchmark):
    oil = benchmark(adios.dirjs.get_full_oil_from_id, 'NO00108')
    assert oil.name == 'AASGARD A 2003'
    assert oil.valid()

def test_extra_oils_no():
    oil = adios.dirjs.find_full_oil_from_name('VEGA CONDENSATE 2015')
    print(oil)

def test_extra_oils_generic():
    oil = adios.dirjs.find_full_oil_from_name('GENERIC HEAVY FUEL OIL')
    print(oil)
