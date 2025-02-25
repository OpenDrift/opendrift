import copy
import requests
import json
import tarfile
import lzma
from importlib.resources import files
from adios_db.models.oil.oil import Oil
from adios_db.computation import gnome_oil

# Oils marked as gnome_suitable, but can not be converted to ADIOS object
blacklist = ['AD02384']

generic_oils = {
        'GENERIC LIGHT CRUDE': {'adios_id': 'GN00006', 'opendrift_id': 'AD04000'},
        'GENERIC MEDIUM CRUDE': {'adios_id': 'GN00007', 'opendrift_id': 'AD04001'},
        'GENERIC HEAVY CRUDE': {'adios_id': 'GN00004', 'opendrift_id': 'AD04002'},
        'GENERIC GASOLINE': {'adios_id': 'GN00003', 'opendrift_id': 'AD04003'},
        'GENERIC FUEL OIL No.2': {'adios_id': 'AD01427', 'opendrift_id': 'AD04006'},
        'GENERIC DIESEL': {'adios_id': 'GN00002', 'opendrift_id': 'AD04007'},
        'GENERIC HOME HEATING OIL': {'adios_id': 'AD02139', 'opendrift_id': 'AD04008'},
        'GENERIC INTERMEDIATE FUEL OIL 180': {'adios_id': 'AD01676', 'opendrift_id': 'AD04009'},
        'GENERIC INTERMEDIATE FUEL OIL 300': {'adios_id': 'AD02428', 'opendrift_id': 'AD04010'},
        'GENERIC FUEL OIL No. 6': {'adios_id': 'AD02431', 'opendrift_id': 'AD04011'},
        'GENERIC BUNKER C': {'adios_id': 'AD02184', 'opendrift_id': 'AD04012'},
        'GENERIC HEAVY FUEL OIL': {'adios_id': 'AD02052', 'opendrift_id': 'AD04013'},
        }

def download():
    #url = 'https://github.com/OpenDrift/noaa-oil-data/archive/refs/heads/production.tar.gz'
    url = 'https://github.com/OpenDrift/noaa-oil-data/archive/refs/heads/new_oils.tar.gz'
    print(f'Downloading all oils from {url}')
    oils = []
    response = requests.get(url, stream=True)
    tar = tarfile.open(fileobj=response.raw, mode='r|gz')
    for file in tar:
        if file.isfile() and file.name.endswith('.json'):
            f = tar.extractfile(file)  
            o = json.load(f.raw)
            AO = Oil.from_py_json(o)
            AO.validate()
            if not AO.metadata.gnome_suitable:
                print(f'Discarding {AO.metadata.name}, not GNOME suitable')
                continue
            else:
                oils.append(o)
    # Add mapping to generic oils
    print('Adding generic oils')
    for generic_name, oil_ids in generic_oils.items():
        adios_id = oil_ids['adios_id']
        opendrift_id = oil_ids['opendrift_id']
        found = False
        for o in oils.copy():
            o = copy.deepcopy(o)
            if o['oil_id'] == adios_id:
                found = True
                break
        if found is False:
            raise ValueError(f'Did not find {adios_id} as {generic_name}')
        else:
            print(f'Mapping {generic_name} [{opendrift_id}] to adios {adios_id}')
        o['metadata']['name'] = generic_name
        #o['metadata']['source_id'] = opendrift_id
        o['oil_id'] = opendrift_id
        oils.append(o)  # Adding modified generic oil

    print(f'Downloaded {len(oils)} oils, saving to oils.xz')
    with lzma.open('oils.xz', 'wt') as c:
        json.dump(oils, c)

def list_oils():
    oil_file = files('opendrift.models.openoil.adios').joinpath('oils.xz')
    with lzma.open(oil_file, 'rt') as archive:
        oils = json.load(archive)

    print(oils)
    print(len(oils))
    for o in oils:
        print(o['oil_id'], o['metadata']['name'])

if __name__ == '__main__':
    download()
    #list_oils()
