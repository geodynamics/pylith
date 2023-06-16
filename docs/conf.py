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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'PyLith'
copyright = '2010-2023, University of California, Davis'
author = 'Brad T. Aagaard, Matthew G. Knepley, Charles A. Williams'

# The full version, including alpha/beta/rc tags
release = '3.0.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
]
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "substitution",
    "dollarmath",
    "amsmath",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = "_static/images/cig_short_pylith.png"
html_theme = 'pydata_sphinx_theme'
html_context = {
    "default_mode": "light",
    "github_user": "geodynamics",
    "github_repo": "pylith",
    "github_version": "main",
    "doc_path": "docs",
}
html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 3,
    "show_toc_level": 3,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/geodynamics/pylith",
            "icon": "fab fa-github-square",
        },
    ],
    "use_edit_page_button": True,
    "navbar_start": ["navbar-logo"],
    "footer_end": ["last-updated"],
    "primary_sidebar_end": ["sidebar-cig"],
}

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "alpha"
bibtex_reference_style = "author_year"
numfig = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
html_last_updated_fmt = ""
