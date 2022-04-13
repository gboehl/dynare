function [ast] = getEquationsByTags(ast, tagname, tagvalue)
%function [ast] = getEquationsByTags(ast, tagname, tagvalue)
% Return the ast structure with the matching tags
%
% INPUTS
%   ast       [cell array]    JSON representation of model block
%   tagname   [string]        The name of the tag whos values are to
%                             be selected
%   tagvalue  [string]        The values to be selected for the
%                             provided tagname
%
% OUTPUTS
%   ast       [cell array]    JSON representation of model block,
%                             with equations removed that don't match
%                             eqtags
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2017-2021 Dynare Team
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

if nargin ~= 3
    error('Incorrect number of arguments passed to function')
end

if isempty(ast) || ~iscell(ast)
    error('the first argument must be a cell array of structs');
end

if ~ischar(tagname)
    error('Tag name must be a string');
end

if ~ischar(tagvalue) && ~iscell(tagvalue)
    error('Tag value must be a string or a cell string array');
end

if ischar(tagvalue)
    tagvalue = {tagvalue};
end

idx2keep = [];
for i = 1:length(tagvalue)
    found = false;
    for j=1:length(ast)
        assert(isstruct(ast{j}), 'Every entry in the ast must be a struct');
        if isfield(ast{j}, 'tags') && ...
                isfield(ast{j}.tags, tagname) && ...
                strcmp(ast{j}.tags.(tagname), tagvalue{i})
            idx2keep = [idx2keep; j];
            found = true;
            break
        end
    end
    if found == false
        warning(['getEquationsByTags: no equation tag found by the name of ''' tagvalue{i} ''''])
    end
end
assert(~isempty(idx2keep), 'getEquationsByTags: no equations selected');
ast = ast(unique(idx2keep, 'stable'));
end
