# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "valladopy"
copyright = "2025, D Vallado, @samrolander"
author = "D Vallado, @samrolander"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "autoapi.extension",
    "myst_nb",  # For Jupyter notebooks
    "sphinx.ext.mathjax",  # For math rendering
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "nbsphinx",  # For supporting .ipynb files (optional)
    "sphinx_collections",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

collections = {
    "notebooks": {
        "driver": "copy_folder",
        "source": "../notebooks",
        "target": "notebooks/",
        "ignore": ["*.py", "*.sh", "*.txt"],
    }
}

autoapi_dirs = ["../src"]
autoapi_root = "API Reference"
autoapi_add_toctree_entry = True
suppress_warnings = ["autoapi"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
    "vcs_pageview_mode": "",
    "flyout_display": "hidden",
    "version_selector": True,
    "language_selector": True,
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}
