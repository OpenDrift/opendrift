from setuptools import setup, find_packages 

config = {
    'name': 'OpenDrift',
    'description': 'OpenDrift - framework for ocean trajectory modeling',
    'author': 'Knut-Frode Dagestad / MET Norway',
    'url': 'https://github.com/knutfrode/opendrift',
    'download_url': 'https://github.com/knutfrode/opendrift',
    'version': '1.0.0',
    'license': 'GPLv2',
    'install_requires': [
        'configobj',
        'numpy',
        'scipy',
        'matplotlib',
        'netCDF4',
        'basemap'
    ],
    'packages': find_packages(),
    'include_package_data': True,
    'test_suite': 'tests',
    #'use_scm_version': True,
    'setup_requires': ['setuptools_scm'],
    'scripts': ['opendrift/scripts/hodograph.py',
                'opendrift/scripts/readerinfo.py',
                'opendrift/scripts/opendrift_plot.py',
                'opendrift/scripts/opendrift_animate.py']
}

setup(**config)
