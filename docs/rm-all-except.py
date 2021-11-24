#! /usr/bin/env python
#
# Remove all empty example files, outside current shard, when building
# sharded-gallery.
#
# Author: Gaute Hope (gauteh@met.no) (c) 2021-09-23

import os
import sys
from pathlib import Path
from glob import glob

with open(sys.argv[1], 'r') as fd:
    files = [Path(f.strip()).stem for f in fd.readlines()]

gallery = Path(__file__).parent / 'source/gallery'
print(f"{files=}")
print(f"{gallery=}")

def check(p: str):
    return not any(p.startswith(f) for f in files)

print('removing example files..')
ex = [Path(e).name for e in glob(str(gallery) + '/example_*')]
ex = list(filter(check, ex))

for e in ex:
    f = gallery / Path(e)
    os.unlink(f)

print('removing thumbnails..')
tex = [Path(e).name[9:] for e in glob(str(gallery) + '/images/thumb/*')]
tex = list(filter(check, tex))

for e in tex:
    f = gallery / Path('images/thumb') / Path('sphx_glr_' + e)
    os.unlink(f)


