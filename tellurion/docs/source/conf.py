# Configuration file for the Sphinx documentation builder.
# pip install sphinx sphinx-rtd-theme sphinx-astropy sphinx-automodapi \
#             sphinx-changelog sphinx-design sphinxcontrib-globalsubs \
#             nbsphinx nbsphinx-link matplotlib ipython pydata-sphinx-theme

import os
import sys
from sphinx_astropy.conf.v2 import *

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
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx_automodapi.automodapi',  # AstroPy's API doc generator
    'sphinx_automodapi.smart_resolver',
    'sphinx_changelog',
    'sphinx_design',
    'sphinxcontrib.globalsubs',
    'nbsphinx',
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
    'Quantity': ': class:`~astropy.units.Quantity`',
    # ... add your project-specific aliases
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
}

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
html_theme = 'pydata_sphinx_theme' # or 'sphinx_rtd_theme', 'alabaster', 'sphinx_book_theme'
# html_theme = 'sphinx_rtd_theme'  # Use Read the Docs theme
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

# This helps with monkey-patched methods
autodoc_mock_imports = []

# Configure automodapi
automodapi_toctreedirnm = 'api'
automodsumm_writereprocessed = False  # Don't write separate reprocessed files
automodapi_writereprocessed = False   # Don't write separate API files
numpydoc_show_class_members = False
autosummary_generate = False  # Don't auto-generate, let automodapi handle it

# Not nitpicky about references
nitpicky = False
