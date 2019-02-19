# -*- coding: utf-8 -*-

# Copyright (C) 2018-2019 Dynare Team
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

sys.path.insert(0, os.path.abspath('../../utils'))

extensions = ['sphinx.ext.autodoc',
			  'sphinx.ext.mathjax']

source_suffix = '.rst'

templates_path = ['_templates']

html_static_path = ['_static']

mathjax_path = 'mathjax/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

master_doc = 'index'

project = u'Dynare'
copyright = u'2019, Dynare Team'
author = u'Dynare Team'

add_function_parentheses = False

version = u'4.6-unstable'
release = u'4.6-unstable'

language = 'en'

exclude_patterns = []

highlight_language = 'dynare'

todo_include_todos = False

html_theme = 'alabaster'

html_sidebars = {
    "**": [
        "about.html",
        "searchbox.html",
        "navigation.html",
    ]
}

html_theme_options = {
    'logo': 'dlogo.svg',
    'logo_name': False,
    'fixed_sidebar': True,
    'page_width': '100%',
}

htmlhelp_basename = 'Dynaremanual'

latex_elements = {
    'sphinxsetup': 'VerbatimBorderColor={rgb}{1,1,1},VerbatimColor={RGB}{240,240,240}, \
                    warningBorderColor={RGB}{255,50,50},OuterLinkColor={RGB}{34,139,34}, \
                    InnerLinkColor={RGB}{51,51,255},TitleColor={RGB}{51,51,255}',
    'papersize': 'a4paper',
}

latex_documents = [
    (master_doc, 'dynare.tex', u'Dynare Reference Manual',
     u'Dynare team', 'manual'),
]

man_pages = [
    (master_doc, 'dynare', u'Dynare Reference Manual',
     [author], 1)
]

def setup(app):
    from dynare_dom import DynareDomain
    from dynare_lex import DynareLexer
    app.add_lexer("dynare", DynareLexer())
    app.add_domain(DynareDomain)
