# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import warnings
sys.path.insert(0, os.path.abspath('../..'))

sys.setrecursionlimit(3000)


# -- Project information -----------------------------------------------------

project = 'OpenDrift'
copyright = '2020, Knut-Frode Dagestad (knutfd@met.no) and Gaute Hope (gauteh@met.no).'
author = 'Knut-Frode Dagestad and Gaute Hope.'
master_doc = 'index'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

autoapi_type = 'python'
autoapi_dirs = [ '../../opendrift' ]
autoapi_keep_files = True  # set to True when debugging autoapi generated files
autoapi_python_class_content = 'both'
autodoc_typehints = 'description'

sphinx_gallery_conf = {
     'examples_dirs': '../../examples/',   # path to your example scripts
     'gallery_dirs': './gallery',  # path to where to save gallery generated output,
     'filename_pattern': '/example_(?!long_)',
     'backreferences_dir': None,
     'capture_repr': ('_repr_html_', '__repr__'),
     'abort_on_example_error': False,
     'thumbnail_size': (300, 300),
     'junit': '../../test-results/sphinx-gallery/junit.xml',
}

# Remove matplotlib agg warnings from generated doc when using plt.show
warnings.filterwarnings("ignore", category=UserWarning,
                        message='Matplotlib is currently using agg, which is a'
                                ' non-GUI backend, so cannot show the figure.')

extensions = [
        "sphinx_rtd_theme",
        "autoapi.extension",
        "sphinx.ext.autodoc",
        "sphinx.ext.viewcode",
        "sphinx.ext.mathjax",
        "sphinx_gallery.gen_gallery",
        "matplotlib.sphinxext.plot_directive",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = '../opendrift_logo.png'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
        'theme_overrides.css'
        ]

