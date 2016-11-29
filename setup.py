from setuptools import setup, find_packages 
import unittest

def opendrift_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite

config = {
    'name': 'OpenDrift',
    'description': 'OpenDrift - framework for ocean trajectory modeling',
    'author': 'Knut-Frode Dagestad / MET Norway',
    'url': 'https://github.com/OpenDrift/opendrift',
    'download_url': 'https://github.com/OpenDrift/opendrift',
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
    'test_suite': 'setup.opendrift_tests',
    #'use_scm_version': True,
    'setup_requires': ['setuptools_scm'],
    'scripts': ['opendrift/scripts/hodograph.py',
                'opendrift/scripts/readerinfo.py',
                'opendrift/scripts/opendrift_plot.py',
                'opendrift/scripts/opendrift_animate.py']
}

setup(**config)
