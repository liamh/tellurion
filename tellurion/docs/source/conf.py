# Configuration file for the Sphinx documentation builder.
# pip install sphinx sphinx-rtd-theme sphinx-astropy sphinx-automodapi \
#             sphinx-changelog sphinx-design sphinxcontrib-globalsubs \
#             myst-nb jupytext matplotlib ipython pydata-sphinx-theme

import os
import sys
from sphinx_astropy.conf.v2 import *

# Compute project root relative to this conf.py file (docs/source/ -> project_root)
project_root = os.path.abspath("../..")

if project_root not in sys.path:
    sys.path.insert(0, project_root)

# -- Project information -----------------------------------------------------
project = 'Tellurion'
copyright = '2026, Liam M. Healy'
author = 'Liam M. Healy'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',        # Auto-generate docs from docstrings
    'sphinx.ext.napoleon',       # Support for numpy-style docstrings
    'sphinx.ext.intersphinx',    # Link to other project docs (like AstroPy)
    'sphinx.ext.viewcode',       # Add links to source code
    'sphinx.ext.mathjax',        # Render math equations
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx_automodapi.automodapi',  # AstroPy's API doc generator
    'sphinx_automodapi.smart_resolver',
    'sphinx_changelog',
    'sphinx_design',
    'sphinxcontrib.globalsubs',
    'myst_nb',                   # Replaced nbsphinx
    'matplotlib.sphinxext.plot_directive',
    'IPython.sphinxext.ipython_console_highlighting',  # For better syntax highlighting
]

suppress_warnings = [
    'app.add_directive',
    'app.add_node',
    'ref.python',
    'autosummary.import_cycle',
    'automodapi',
    'toc.not_included',
    'ref.ref',
    'autosummary',
    'autodoc',  # Suppress "don't know which module" warnings
    'docutils',  # Suppress the "Explicit markup ends without blank line" warning
]

plot_rcparams = {
    'axes.labelsize': 'large',
    'figure.figsize': (6, 6),
    'savefig.bbox': 'tight',
}
plot_apply_rcparams = True
plot_html_show_source_link = False
plot_formats = ['png', 'svg', 'pdf']

numpydoc_xref_param_type = True
numpydoc_xref_aliases = {
    'Quantity': ':class:`~astropy.units.Quantity`',
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
}

# -- myst-nb & Jupytext Configuration ----------------------------------------
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
    '.md': 'myst-nb',
    '.py': 'myst-nb',
}

myst_heading_anchors = 3  # Automatically creates anchors for h1, h2, and h3 headers

jupytext_custom_formats = {
    ".py": "jupytext.reads(RST_TEXT, format_name='percent')",
}

# Ensure MyST-NB recognizes jupytext percent files
nb_custom_formats = {
    ".py": ["jupytext.reads", {"fmt": "py:percent"}]
}

nb_kernel_args = ["--transport=ipc"]
nb_execution_mode = 'auto'
nb_execution_timeout = 180

# Napoleon settings for numpy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_ivar = True

# Automodapi settings
numpydoc_show_class_members = False

templates_path = ['_templates']
exclude_patterns = [
    "tutorials/cartprop2.py",
]

# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_theme_options = dict(globals().get('html_theme_options', {}))
html_theme_options['secondary_sidebar_items'] = []
html_theme_options['show_toc_level'] = 4
html_sidebars = {
    '**': ['sidebar-nav-bs', 'page-toc'],
}

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

autodoc_mock_imports = []

# Configure automodapi
automodapi_toctreedirnm = 'api'
automodsumm_writereprocessed = False
automodapi_writereprocessed = False
autosummary_generate = False

nitpicky = False
