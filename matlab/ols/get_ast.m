function [ast, jsonmodel] = get_ast(eqtags)

% Returns the Abstract Syntax Tree the JSON output of preprocessor for the
% given equation tags. Equations are ordered in eqtag order
%
% INPUTS
% - eqtags     [cellstr]    names of equation tags for which to get info.
%                           If empty, get all equations
%
% OUTPUTS
% - ast        [cell array]  JSON representation of the abstract syntax tree
% - jsonmodel  [cell array]  JSON representation of the equations.
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

global M_

if ~isempty(eqtags) && ~(ischar(eqtags) || iscell(eqtags))
    error('eqtags not passed correctly')
end

jsonfile = sprintf('%s%smodel%sjson%smodfile-original.json', M_.fname, filesep(), filesep(), filesep());

if ~exist(jsonfile, 'file')
    error('Could not find %s! Please use the json=compute option (see the Dynare invocation section in the reference manual).', jsonfile)
end

tmp = loadjson_(jsonfile);
ast = tmp.abstract_syntax_tree;

if nargout>1
    jsonmodel = tmp.model;
end

if ~isempty(eqtags)
    ast = getEquationsByTags(ast, 'name', eqtags);
    if nargout>1
        jsonmodel = getEquationsByTags(jsonmodel, 'name', eqtags);
    end
end