# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "valladopy"
copyright = "2026, D. A. Vallado and Contributors"
author = "D. A. Vallado, maintained by S. Rolander"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "autoapi.extension",
    "sphinx.ext.mathjax",  # For math rendering
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "nbsphinx",  # For supporting .ipynb files
    "sphinx_collections",
    "sphinx_copybutton",
]
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

collections = {
    "notebooks": {
        "driver": "copy_folder",
        "source": "../notebooks",
        "target": "notebooks/",
        "ignore": ["*.py", "*.sh"],
    }
}

nbsphinx_kernel_name = "vallado-env"


autoapi_dirs = ["../src"]
autoapi_root = "API Reference"
autoapi_add_toctree_entry = True
suppress_warnings = ["autoapi"]

html_static_path = ["_static"]

html_theme = "sphinx_book_theme"
html_title = "valladopy"
html_copy_source = False
html_show_sourcelink = False
# furo
html_theme_options = {
    "top_of_page_buttons": [],
}
# sphinx book
html_theme_options = {
    "repository_url": "https://github.com/celestrak/fundamentals-of-astrodynamics",
    "use_repository_button": True,
}
