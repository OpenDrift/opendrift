import requests
import json
import tarfile
import lzma
from adios_db.models.oil.oil import Oil
from adios_db.computation import gnome_oil

# Oils marked as gnome_suitable, but can not be converted to ADIOS object
blacklist = ['AD02384']

def download():
    url = 'https://github.com/OpenDrift/noaa-oil-data/archive/refs/heads/production.tar.gz'
    print(f'Downloading all oils from {url}')
    oils = []
    response = requests.get(url, stream=True)
    tar = tarfile.open(fileobj=response.raw, mode='r|gz')
    for file in tar:
        if file.isfile() and file.name.endswith('.json'):
            gnome_suitable = True
            for b in blacklist:
                if b in file.name:
                    print(f'Discarding {b}, from blacklist (not GNOME suitable)')
                    gnome_suitable = False
            if gnome_suitable is False:
                continue
            f = tar.extractfile(file)  
            o = json.load(f.raw)
            AO = Oil.from_py_json(o)
            #if 'norw' in AO.metadata.location.lower() and not AO.metadata.gnome_suitable:
            if not AO.metadata.gnome_suitable:
                print(f'Discarding {AO.metadata.name}, not GNOME suitable')
                continue
            else:
                oils.append(o)
    print(f'Downloaded {len(oils)} oils, saving to oils.xz')
    with lzma.open('oils.xz', 'wt') as c:
        json.dump(oils, c)

def list_oils():
    from importlib import resources
    import lzma
    with resources.open_binary('opendrift.models.openoil.adios',
        'oils.xz') as archive:
        with lzma.open(archive, 'rt') as c:
            oils = json.load(c)
    print(oils)
    print(len(oils))
    for o in oils:
        print(o['oil_id'], o['metadata']['name'])

if __name__ == '__main__':
    download()
    #list_oils()
