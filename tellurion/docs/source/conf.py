# Configuration file for the Sphinx documentation builder.

import os
import sys

# Add your package to the path so Sphinx can import it
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'Tellurion'
copyright = '2025, Liam M. Healy'
author = 'Liam Healy'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',        # Auto-generate docs from docstrings
    'sphinx.ext.napoleon',       # Support for numpy-style docstrings
    'sphinx.ext.intersphinx',    # Link to other project docs (like AstroPy)
    'sphinx.ext.viewcode',       # Add links to source code
    'sphinx.ext.mathjax',        # Render math equations
    'sphinx_automodapi.automodapi',  # AstroPy's API doc generator
    'sphinx_automodapi.smart_resolver',
    'nbsphinx',
    'nbsphinx_link',
    'IPython.sphinxext.ipython_console_highlighting',  # For better syntax highlighting
]

# Suppress the harmless cache warning
suppress_warnings = ['config.cache']

# Don't execute notebooks during build (use pre-executed outputs)
nbsphinx_execute = 'never'

# Timeout for notebook execution (if you do execute)
nbsphinx_timeout = 180

# Napoleon settings for numpy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

# Intersphinx mapping to link to external docs
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# Automodapi settings
numpydoc_show_class_members = False

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# html_theme = 'sphinx_astropy' # or 'sphinx_rtd_theme', 'alabaster', 'sphinx_book_theme'
html_theme = 'sphinx_rtd_theme'  # Use Read the Docs theme
html_static_path = ['_static']
