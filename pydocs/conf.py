# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GFEAT'
copyright = '2025, Gabriel Valles Valverde'
author = 'Gabriel Valles Valverde'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",           
    "sphinx.ext.viewcode",           
    "sphinx_autodoc_typehints",  
    'sphinx.ext.inheritance_diagram',
    'nbsphinx',    
    'enum_tools.autoenum',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

nbsphinx_thumbnails = {
    'notebooks/gravity-field-data/aod1b': '../../notebooks/gravity-field-data/aod1b-2025-01.gif',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']


def autodoc_skip_member(app, what, name, obj, skip, options):
    if name == 'pybind11_object':
        return True
    return None


def setup(app):
    app.add_css_file('custom.css')
    app.connect('autodoc-skip-member', autodoc_skip_member)

autoclass_content = "both"
