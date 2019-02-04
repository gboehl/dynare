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

"""
    sphinx.domains.dynare
    ~~~~~~~~~~~~~~~~~~~~~

    The Dynare domain.

    Loosely based on the JavaScript domain from the Sphinx team and the CoffeScript domain
    by Stephen Sugden (available at sphinx-contrib)
"""

import re

from docutils import nodes
from docutils.parsers.rst import Directive, directives

from sphinx import addnodes
from sphinx.domains import Domain, ObjType
from sphinx.locale import l_, _
from sphinx.directives import ObjectDescription
from sphinx.roles import XRefRole
from sphinx.util.nodes import make_refnode
from sphinx.util.docfields import Field, GroupedField, TypedField

############### Dynare Object types #######################

class DynObject(ObjectDescription):
    has_arguments  = True
    display_prefix = None
    allow_nesting  = False          

    def handle_signature(self, sig, signode):
        sig = sig.strip()

        # Some variable/parameter declarations combine spaces and parenthesis
        # in a strange way, so they are treated separately
        if sig.startswith(('var V','varexo V','varexo_det V','parameters P','model_comparison F')):
            member  = sig[:sig.index(' ')]
            arglist = sig[sig.index(' '):]
        # General cases
        elif '(' in sig:
            member  = sig[:sig.index('(')]
            arglist = sig[sig.index('('):]
            arglist = arglist.strip()
        elif ' ' in sig:
            member  = sig[:sig.index(' ')]
            arglist = sig[sig.index(' '):]
        else:
            member = sig
            arglist = None

        prefix = self.env.ref_context.get('dynare:object', None)
        name = member.strip()
        fullname = name

        signode['object'] = prefix
        signode['fullname'] = fullname

        if self.display_prefix:
            signode += addnodes.desc_annotation(self.display_prefix, self.display_prefix)
        
        signode += addnodes.desc_name(name, name)
        
        if self.has_arguments:
            if not arglist:
                signode += addnodes.desc_parameterlist()
            else:
                signode += addnodes.desc_addname(arglist,arglist)
        return fullname, prefix

    def add_target_and_index(self, name_obj, sig, signode):
        fullname = name_obj[0]
        if fullname not in self.state.document.ids:
            signode['names'].append(fullname)
            signode['ids'].append(fullname)
            signode['first'] = not self.names
            self.state.document.note_explicit_target(signode)
            objects = self.env.domaindata['dynare']['objects']
            
            if fullname in objects:
                self.state_machine.reporter.warning(
                    'duplicate object description of %s, ' % fullname +
                    'other instance in ' +
                    self.env.doc2path(objects[fullname][0]),line=self.lineno)
            objects[fullname] = (self.env.docname, self.objtype)

        indextext = self.get_index_text(fullname,name_obj)
        if indextext:
            self.indexnode['entries'].append(('single', indextext, fullname,'', None))

    def get_index_text(self, objectname, name_obj):
        name, obj = name_obj

        regexp = re.compile(r'\s*=\s*')
        if bool(regexp.search(name)):
            aux, name = re.compile(r'\s*=\s*').split(name)

        if self.objtype == 'function':
            return _('%s (function)') % name
        elif self.objtype == 'class':
            return _('%s (class)') % name
        elif self.objtype == 'datesmethod':
            return _('%s (dates method)') % name
        elif self.objtype == 'dseriesmethod':
            return _('%s (dseries method)') % name
        elif self.objtype == 'reportingmethod':
            return _('%s (reporting method)') % name
        elif self.objtype == 'matcomm':
            return _('%s (MATLAB command)') % name
        elif self.objtype == 'command':
            return _('%s (command)') % name
        elif self.objtype == 'block':
            return _('%s (block)') % name
        elif self.objtype == 'confblock':
            name = name[1:-1]
            return _('%s (config block)') % name
        elif self.objtype == 'macrodir':
            name = name[2:]
            return _('%s (macro directive)') % name

class DynCallable(DynObject):
    has_arguments = True

    doc_field_types = [
        TypedField('arguments', label=l_('Arguments'),
                   names=('argument', 'arg', 'parameter', 'param'),
                   typerolename='func', typenames=('paramtype', 'type')),
        Field('returnvalue', label=l_('Returns'), has_arg=False,
              names=('returns', 'return')),
        Field('returntype', label=l_('Return type'), has_arg=False,
              names=('rtype',)),
        Field('example', label=l_('Example'), has_arg=False,
              names=('ex',)),
    ]

class DynClass(DynObject):
    has_arguments = False
    display_prefix = 'Dynare class: '
    allow_nesting = True

    doc_field_types = [
        TypedField('members', label=l_('Members'),
                   names=('argument', 'arg', ),
                   typerolename='func', typenames=('type', )),
        Field('example', label=l_('Example'), has_arg=False,
              names=('ex',)),
    ]

class DynFunction(DynCallable):
    display_prefix = 'Function: '
    allow_nesting = True

class DatesMethod(DynCallable):
    display_prefix = 'Method: '
    allow_nesting = True

class DseriesMethod(DynCallable):
    display_prefix = 'Method: '
    allow_nesting = True

class ReportingMethod(DynCallable):
    display_prefix = 'Method: '
    allow_nesting = True

class MatComm(DynCallable):
    display_prefix = 'MATLAB/Octave command: '
    allow_nesting = False

class DynComm(DynCallable):
    display_prefix = 'Command: '
    allow_nesting = False

class DynBlock(DynCallable):
    display_prefix = 'Block: '
    allow_nesting = False

class DynConfBlock(DynCallable):
    display_prefix = 'Configuration block: '
    has_arguments = False
    allow_nesting = False

class DynMacroDir(DynCallable):
    display_prefix = 'Macro directive: '
    allow_nesting = False

class Constructor(DynCallable):
    display_prefix = 'Constructor: '
    allow_nesting = False

class DynSimpleObject(ObjectDescription):
    has_arguments = False
    allow_nesting = False

    def handle_signature(self, sig, signode):
        sig = sig.strip()
        member = sig
        arglist = None
        prefix = self.env.ref_context.get('dynare:object', None)
        name = member
        fullname = name

        signode['object'] = prefix
        signode['fullname'] = fullname

        if self.display_prefix:
            signode += addnodes.desc_annotation(self.display_prefix, self.display_prefix)
        
        signode += addnodes.desc_name(name, name)
        return fullname, prefix    

    def add_target_and_index(self, name_obj, sig, signode):
        fullname = name_obj[0]
        if fullname not in self.state.document.ids:
            signode['names'].append(fullname)
            signode['ids'].append(fullname)
            signode['first'] = not self.names
            self.state.document.note_explicit_target(signode)
            objects = self.env.domaindata['dynare']['objects']
            if fullname in objects:
                self.state_machine.reporter.warning(
                    'duplicate object description of %s, ' % fullname +
                    'other instance in ' +
                    self.env.doc2path(objects[fullname][0]), line=self.lineno)
            objects[fullname] = self.env.docname, self.objtype

        indextext = self.get_index_text(fullname,name_obj)
        if indextext:
            self.indexnode['entries'].append(('single', indextext, fullname,'', None))

    def get_index_text(self, objectname, name_obj):
        name, obj = name_obj

        if self.objtype == 'construct':
            name, rest = name.split(' ',1)
            return _('%s (constructor)') % name
        elif self.objtype == 'matvar':
            return _('%s (MATLAB variable)') % name
        elif self.objtype == 'specvar':
            return _('%s (special variable)') % name
        elif self.objtype == 'operator':
            endsig = name.find(' ')
            name = name[0:endsig]
            return _('%s (operator)') % name
        elif self.objtype == 'constant':
            return _('%s (constant)') % name

class MatlabVar(DynSimpleObject):
    display_prefix = 'MATLAB/Octave variable: '
    allow_nesting = False

class SpecialVar(MatlabVar):
    display_prefix = 'Special variable: '

class Operator(MatlabVar):
    display_prefix = 'Operator: '

class Constant(MatlabVar):
    display_prefix = 'Constant: '

class Option(MatlabVar):
    display_prefix = None

############## Cross-referencing ####################

class DynareXRefRole(XRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['dynare:object'] = env.ref_context.get('dynare:object')
        return title, target

############### Dynare domain #######################

class DynDomain(Domain):
    name = 'dynare'
    label = 'Dynare'
    object_types = {
    'function':      ObjType(l_('function'),         'func'),
    'datesmethod':   ObjType(l_('method'),           'datmeth'),
    'dseriesmethod': ObjType(l_('method'),           'dsermeth'),
    'reportingmethod': 	 ObjType(l_('method'),       'repmeth'),
    'matcomm':   	 ObjType(l_('matlab command'),   'mcomm'),
    'command':   	 ObjType(l_('command'),          'comm'),
    'class':     	 ObjType(l_('class'),            'class'),
    'block':     	 ObjType(l_('block'),            'bck'),
    'confblock': 	 ObjType(l_('config block'),     'cbck'),
    'macrodir':  	 ObjType(l_('macro directive'),  'mdir'),
    'construct': 	 ObjType(l_('constructor'),      'cstr'),
    'matvar':    	 ObjType(l_('matlab variable'),  'mvar'),
    'specvar':   	 ObjType(l_('special variable'), 'svar'),
    'operator':  	 ObjType(l_('operator'),         'op'),
    'constant':  	 ObjType(l_('constant'),         'const'),
    'option':    	 ObjType(l_('option'),           'opt'),
    }
    directives = {
    'function':      DynFunction,
    'datesmethod':   DatesMethod,
    'dseriesmethod': DseriesMethod,
    'reportingmethod': 	 ReportingMethod,
    'matcomm':   	 MatComm,
    'command':   	 DynComm,
    'class':     	 DynClass,
    'block':     	 DynBlock,
    'confblock': 	 DynConfBlock,
    'macrodir':  	 DynMacroDir,
    'construct': 	 Constructor,   
    'matvar':    	 MatlabVar,
    'specvar':   	 SpecialVar,  
    'operator':  	 Operator, 
    'constant':  	 Constant,
    'option':    	 Option,
    }
    roles = {
    'func':     DynareXRefRole(),
    'datmeth':  DynareXRefRole(),
    'dsermeth': DynareXRefRole(),
    'repmeth':  DynareXRefRole(),
    'mcomm':    DynareXRefRole(),
    'comm':     DynareXRefRole(),
    'class':    DynareXRefRole(),
    'bck':      DynareXRefRole(),    
    'cbck':     DynareXRefRole(),
    'mdir':     DynareXRefRole(),
    'cstr':     DynareXRefRole(),       
    'mvar':     DynareXRefRole(),   
    'svar':     DynareXRefRole(),    
    'op':       DynareXRefRole(),     
    'const':    DynareXRefRole(),
    'opt':      DynareXRefRole(),  
    }
    initial_data = {
        'objects': {},
    } 

    def clear_doc(self, docname):
        for fullname, (fn, _l) in list(self.data['objects'].items()):
            if fn == docname:
                del self.data['objects'][fullname]

    def find_obj(self, env, obj, name, typ, searchorder=0):
        objects = self.data['objects']
        newname = name # None
        return newname, objects.get(newname)

    def merge_domaindata(self, docnames, otherdata):
        for fullname, (fn, objtype) in otherdata['objects'].items():
            if fn in docnames:
                self.data['objects'][fullname] = (fn, objtype)

    def resolve_xref(self, env, fromdocname, builder, typ, target, node, contnode):
        objectname = node.get('dynare:object')
        searchorder = node.hasattr('refspecific') and 1 or 0
        name, obj = self.find_obj(env, objectname, target, typ, searchorder)
        if not obj:
            return None
        return make_refnode(builder, fromdocname, obj[0], name, contnode, name)

    def get_objects(self):
        for refname, (docname, type) in list(self.data['objects'].items()):
            yield refname, refname, type, docname, refname, 1
