#!/usr/bin/env python
"""
Leeway crosswind fix: before/after maps
========================================

Runs identical Leeway simulations with the Leeway module before the fix of
crosswind direction and jibing (OpenDrift commit ed80c33c) and with the current
code, and plots maps side by side, with trajectories coloured by orientation
(red = right of downwind, blue = left of downwind).

Two cases, both for PIW-4 (person in survival suit, asymmetric crosswind coefficients):

- idealised: constant 10 m/s wind towards north, no current, 24 h, no jibing
- norkyst_arome: NorKyst800 currents and AROME wind from the test data, 48 h

Usage::

    python examples/leeway_crosswind_before_after.py [--old path/to/old/leeway.py] [--outdir leeway_maps]

If --old is not given, the old leeway.py is taken with ``git show`` from the
local repository, or downloaded from GitHub.
"""

import argparse
import importlib.util
import os
import subprocess
import tempfile
import urllib.request
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import opendrift
from opendrift import test_data_folder as tdf
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway as LeewayNew

OLD_COMMIT = 'ed80c33c27d14b29014fde5bfa0aafacfa45a109'  # before the fix
OLD_URL = ('https://raw.githubusercontent.com/OpenDrift/opendrift/'
           f'{OLD_COMMIT}/opendrift/models/leeway.py')
OBJECTPROP = os.path.join(os.path.dirname(opendrift.__file__), 'models', 'OBJECTPROP.DAT')
CMAP = ListedColormap(['tab:red', 'tab:blue'])  # orientation 0: right, 1: left of downwind


def get_old_leeway(path=None):
    """Return the Leeway class from the module before the fix."""
    if path is None:
        path = os.path.join(tempfile.mkdtemp(), 'leeway_old.py')
        repo = os.path.join(os.path.dirname(opendrift.__file__), '..')
        try:
            source = subprocess.run(
                ['git', '-C', repo, 'show', f'{OLD_COMMIT}:opendrift/models/leeway.py'],
                check=True, capture_output=True).stdout
            print(f'Old leeway.py taken from local git, commit {OLD_COMMIT[:8]}')
        except (OSError, subprocess.CalledProcessError):
            source = urllib.request.urlopen(OLD_URL).read()
            print(f'Old leeway.py downloaded from {OLD_URL}')
        with open(path, 'wb') as f:
            f.write(source)
    spec = importlib.util.spec_from_file_location('leeway_old', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Leeway


def idealised(cls):
    o = cls(d=OBJECTPROP, loglevel=50)
    o.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 0, 'y_wind': 10, 'land_binary_mask': 0})
    o.seed_elements(lon=3, lat=59, time=datetime(2020, 1, 1), number=1000,
                    object_type=4, jibe_probability=0)
    o.run(duration=timedelta(hours=24), time_step=900, time_step_output=3600)
    return o


def norkyst_arome(cls):
    o = cls(d=OBJECTPROP, loglevel=50)
    o.add_reader([
        reader_netCDF_CF_generic.Reader(tdf + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc'),
        reader_netCDF_CF_generic.Reader(tdf + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')])
    o.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                    time=datetime(2015, 11, 16, 0), object_type=4)  # default jibing
    o.run(duration=timedelta(hours=48), time_step=900, time_step_output=3600)
    return o


CASES = (
    ('idealised', idealised, .1,
     'PIW-4, constant 10 m/s wind towards north, 24 h, no jibing'),
    ('norkyst_arome', norkyst_arome, .3,
     'PIW-4, NorKyst800 + AROME, 16 Nov 2015, 48 h, jibing 0.04/h'),
)


def main(old=None, outdir='leeway_maps'):
    os.makedirs(outdir, exist_ok=True)
    LeewayOld = get_old_leeway(old)

    for name, func, pad, title in CASES:
        o_old, o_new = func(LeewayOld), func(LeewayNew)

        # Same map extent for before and after
        lons = np.concatenate([o.result.lon.values.ravel() for o in (o_old, o_new)])
        lats = np.concatenate([o.result.lat.values.ravel() for o in (o_old, o_new)])
        corners = [np.nanmin(lons) - pad, np.nanmax(lons) + pad,
                   np.nanmin(lats) - pad / 3, np.nanmax(lats) + pad / 3]

        files = []
        for tag, o in (('before', o_old), ('after', o_new)):
            filename = os.path.join(outdir, f'{name}_{tag}.png')
            o.plot(fast=True, linecolor='orientation', cmap=CMAP, lvmin=-.5, lvmax=1.5,
                   linewidth=.6, corners=corners, colorbar=False, show=False,
                   title=f'{tag.upper()} fix\nred = right of downwind, blue = left of downwind',
                   filename=filename)
            files.append(filename)
            print(f'{name} {tag}: {o.num_elements_active()} active, '
                  f'{o.num_elements_deactivated()} stranded')

        fig, axes = plt.subplots(1, 2, figsize=(16, 9))
        for ax, filename in zip(axes, files):
            ax.imshow(plt.imread(filename))
            ax.axis('off')
        fig.suptitle(title, fontsize=15)
        fig.tight_layout()
        filename = os.path.join(outdir, f'{name}_before_after.png')
        fig.savefig(filename, dpi=110)
        plt.close(fig)
        print(f'Saved {filename}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    parser.add_argument('--old', help='Path to leeway.py before the fix '
                        '(default: from local git or GitHub)')
    parser.add_argument('--outdir', default='leeway_maps', help='Output folder for figures')
    args = parser.parse_args()
    main(args.old, args.outdir)
