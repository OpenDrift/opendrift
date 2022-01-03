from pathlib import Path
import os
import sys
import requests
from opendrift.models.openoil.adios.oil import OpendriftOil

import coloredlogs
import logging
logger = logging.getLogger('harvest_oils')

ADIOS = "https://adios.orr.noaa.gov/api/oils/"

def get_full_oil_from_id(_id):
    logger.debug(f"Fetching full oil: {_id}")
    o = requests.get(f"{ADIOS}/{_id}")
    o.raise_for_status()
    return OpendriftOil(o.json())

def oils(limit=9999, query=''):
    """
    Get all oils.

    Args:

        limit: number of oils to retrieve, None means all available.
        query: search string (name, id, location).


    Returns:

        List of `class:ThinOil`s.
    """

    o = requests.get(ADIOS, {
        'dir': 'asc',
        'limit': limit,
        'page': 0,
        'sort': 'metadata.name',
        'q': query
    })
    o.raise_for_status()
    o = o.json()

    oils = o['data']

    logging.debug(f"Total oils: {o['meta']['total']}, download limit: {limit}")
    logger.info(f"Fetched {len(oils)} oils from ADIOS")

    limit = len(oils) if limit is None else limit
    oils = oils[:min(limit, len(oils))]
    oils = [o['_id'] for o in oils]

    return oils

def download(dst):
    logger.info('downloading all oils..')

    ols = oils()
    logger.info(f'downloaded list of oils: {len(ols)}, fetching full oil..')

    for o in ols:
        logger.info(f"fetching oil: {o}..")
        o = get_full_oil_from_id(o)
        if o.valid():
            f = dst / Path(o.id).with_suffix('.json')
            with open(f, 'w') as fd:
                fd.write(o.json())
        else:
            logger.warning(f'skipping invalid oil: {o.id}')

def make_archive(oildir, file):
    logger.info(f'making archive: {oildir / "*.json"}')

    # Making a big JSON array with oils as a dictionary each.
    oils = []

    import glob
    import json
    import lzma
    for f in glob.glob(str(oildir / '*.json')):
        with open(f, 'r') as fd:
            oils.append(json.load(fd))

    logger.info(f'added {len(oils)} oils..')
    logger.info(f'compressing to {file}..')
    with lzma.open(file, 'wt') as c:
        json.dump(oils, c)


if __name__ == '__main__':
    coloredlogs.install('debug')

    dst = Path('oils')

    if not dst.exists():
        os.makedirs(dst)

    if not '--skip-download' in sys.argv:
        download(dst)

    make_archive(dst, 'oils.xz')

