from setuptools import setup, find_packages 

config = {
    'name': 'opendrift',
    'description': 'opendrift - framework for ocean trajectory modeling',
    'author': 'Knut-Frode Dagestad / MET Norway',
    'url': 'https://github.com/knutfrode/opendrift',
    'download_url': 'https://github.com/knutfrode/opendrift',
    'version': '1.0.0',
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
    #'use_scm_version': True,
    'setup_requires': ['setuptools_scm'],
    'scripts': ['scripts/hodograph.py', 'scripts/readerinfo.py',
                'scripts/opendrift_plot.py', 'scripts/opendrift_animate.py']
}

setup(**config)
