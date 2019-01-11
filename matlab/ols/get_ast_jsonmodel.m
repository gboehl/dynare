function [ast, jsonmodel] = get_ast_jsonmodel(eqtags)
% function [ast, jsonmodel] = get_ast_jsonmodel(eqtags)
% return AST and jsonmodel from json output of preprocessor for the given
% equation tags.
%
% INPUTS
%   eqtags     [cellstr]    names of equation tags for which to get info.
%                           If empty, get all
%
% OUTPUTS
%   ast        [cell array] JSON representation of the abstract syntax tree
%   jsonmodel  [cell array] JSON representation of the model block
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_

if ~isempty(eqtags) && ~(ischar(eqtags) || iscell(eqtags))
    error('get_ast_jsonmodel: eqtags not passed correctly')
end

jsonfile = [M_.fname filesep() 'model' filesep() 'json' filesep() 'modfile-original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end

jsonmodel = loadjson(jsonfile);
ast       = jsonmodel.abstract_syntax_tree;
jsonmodel = jsonmodel.model;
if ~isempty(eqtags)
    ast = getEquationsByTags(ast, 'name', eqtags);
    jsonmodel = getEquationsByTags(jsonmodel, 'name', eqtags);
end
end