function [jsonmodel] = getEquationsByTags(jsonmodel, tagname, tagvalue)
%function [jsonmodel] = getEquationsByTags(jsonmodel, tagname, tagvalue)
% Return the jsonmodel structure with the matching tags
%
% INPUTS
%   jsonmodel       [cell array]    JSON representation of model block
%   tagname         [string]        The name of the tag whos values are to
%                                   be selected
%   tagvalue        [string]        The values to be selected for the
%                                   provided tagname
%
% OUTPUTS
%   jsonmodel       [cell array]    JSON representation of model block,
%                                   with equations removed that don't match
%                                   eqtags
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017-2018 Dynare Team
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

assert(nargin == 3, 'Incorrect number of arguments passed to getEquationsByTags');
assert(iscell(jsonmodel) && ~isempty(jsonmodel), ...
    'the first argument must be a cell array of structs');
assert(ischar(tagname), 'Tag name must be a string');
assert(ischar(tagvalue) || iscell(tagvalue), 'Tag value must be a string or a cell string array');

if ischar(tagvalue)
    tagvalue = {tagvalue};
end

idx2keep = [];
for i=1:length(tagvalue)
    found = false;
    for j=1:length(jsonmodel)
        assert(isstruct(jsonmodel{j}), 'Every entry in jsonmodel must be a struct');
        if isfield(jsonmodel{j}, 'tags') && ...
                isfield(jsonmodel{j}.tags, tagname) && ...
                strcmp(jsonmodel{j}.tags.(tagname), tagvalue{i})
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
jsonmodel = jsonmodel(unique(idx2keep, 'stable'));
end
