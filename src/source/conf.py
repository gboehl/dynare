# -*- coding: utf-8 -*-

# Copyright (C) 2018 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys

sys.path.insert(0, os.path.abspath('../../py/domain/'))

extensions = ['sphinx.ext.autodoc',
			  'sphinx.ext.mathjax']

templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

project = u'Dynare'
copyright = u'2018, Dynare Team'
author = u'Dynare Team'

add_function_parentheses = False

version = u'4.6-unstable'
release = u'4.6-unstable'

language = 'en'

exclude_patterns = []

highlight_language = 'dynare'

todo_include_todos = False

html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

html_logo = 'dynarelogo.png'

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Dynaremanual'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {

    'sphinxsetup': 'VerbatimBorderColor={rgb}{1,1,1},VerbatimColor={RGB}{240,240,240}, \
                    warningBorderColor={RGB}{255,50,50},OuterLinkColor={RGB}{34,139,34}, \
                    InnerLinkColor={RGB}{51,51,255},TitleColor={RGB}{51,51,255}',

    'papersize': 'a4paper',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'dynare.tex', u'Dynare Reference Manual',
     u'Dynare team', 'manual'),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'manual', u'Dynare Reference Manual',
     [author], 1)
]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Dynare', u'Dynare Reference Manual',
     author, 'Test', 'One line description of project.',
     'Miscellaneous'),
]

def setup(sphinx):
	from dynare import DynDomain
	sphinx.add_domain(DynDomain)
