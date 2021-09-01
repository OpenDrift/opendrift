from opendrift.models.openoil import adios
from opendrift.models.openoil.adios.oil import OpendriftOil
from pathlib import Path
import shutil
import sys

import requests
requests.packages.urllib3.disable_warnings()

import coloredlogs
import logging
logger = logging.getLogger('harvest_oils')
coloredlogs.install('debug')

dst = Path('oils')

if not dst.exists():
    raise Exception(f"destination path does not exist: {dst}")

if not '--skip-download' in sys.argv:
    for o in adios.oils(None):
        logger.info(f"fetching oil: {o.id} / {o.name}")
        o = o.make_full()
        if o.valid():
            f = dst / Path(o.id).with_suffix('.json')
            with open(f, 'w') as fd:
                fd.write(o.json())
        else:
            logger.warning(f'skipping invalid oil: {o.id}')

logger.info('making archive')
## Tarfile
# import glob
# import tarfile
# with tarfile.open('oils.tar.xz', 'x:xz') as c:
#     for f in glob.glob('oils/*.json'):
#         logger.debug(f'adding {f} to archive')
#         c.add(f, Path(f).name)

## Big json
oils = []

import glob
import json
import lzma
for f in glob.glob('oils/*.json'):
    with open(f, 'r') as fd:
        oils.append(json.load(fd))

logger.info(f'added {len(oils)} oils..')
with lzma.open('oils.xz', 'xt') as c:
    json.dump(oils, c)

## pickle
# import glob
# import pickle
# import json
# oils = []
# for f in glob.glob('oils/*.json'):
#     with open(f, 'r') as fd:
#         oils.append(OpendriftOil(json.load(fd)))

# logger.info(f'added {len(oils)} oils..')
# with open('oils.db', 'wb') as fd:
#     pickle.dump(oils, fd)
