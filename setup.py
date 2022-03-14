#!/usr/bin/env python

import os
import setuptools

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'opendrift/version.py')).read())

setuptools.setup(
    name        = 'OpenDrift',
    description = 'OpenDrift - a framework for ocean trajectory modeling',
    author      = 'Knut-Frode Dagestad / MET Norway',
    url         = 'https://github.com/OpenDrift/opendrift',
    download_url = 'https://github.com/OpenDrift/opendrift',
    version = __version__,
    license = 'GPLv2',
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'netCDF4',
        'pyproj',
        'cartopy',
    ],
    packages = setuptools.find_packages(),
    include_package_data = True,
    setup_requires = ['setuptools_scm'],
    tests_require = ['pytest'],
    scripts = ['opendrift/scripts/hodograph.py',
               'opendrift/scripts/readerinfo.py',
               'opendrift/scripts/opendrift_plot.py',
               'opendrift/scripts/opendrift_animate.py',
               'opendrift/scripts/opendrift_gui.py',
               'opendrift/scripts/mp4_to_gif.bash']
)

