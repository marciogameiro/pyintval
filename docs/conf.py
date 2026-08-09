"""Sphinx configuration for the pyintval documentation."""

import importlib.metadata

project = "pyintval"
author = "Marcio Gameiro"
copyright = "2026, Marcio Gameiro"

try:
    release = importlib.metadata.version("pyintval")
except importlib.metadata.PackageNotFoundError:
    release = "0.1.0.dev0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True}
napoleon_google_docstring = True

myst_enable_extensions = ["dollarmath", "colon_fence"]

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = f"pyintval {release}"
