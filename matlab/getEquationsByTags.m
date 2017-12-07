function [lhs, rhs, linenum, sample, tagvalue] = getEquationsByTags(jsonmodel, varargin)
%function [lhs, rhs, linenum, sample] = getEquationByTag(jsonmodel, varargin)
% Return the lhs, rhs of an equation and the line it was defined
% on given its tag
%
% INPUTS
%   jsonmodel        [string] JSON representation of model block
%   varargin         [string or cellstring arrays] tagname and tagvalue for
%                                                  eqs to get
%
% OUTPUTS
%   lhs             [cellstring array]     left hand side of eq
%   rhs             [cellstring array]     right hand side of eq
%   linenum         [cellstring array]     eq line in .mod file
%   sample          [cell array of dates]  sample range
%   tagvalue        [cellstring array]     tags associated with equations
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017 Dynare Team
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

assert(nargin == 1 || nargin == 3, 'Incorrect number of arguments passed to getEquationsByTags');

if nargin == 1
    lhs = cell(1, length(jsonmodel));
    rhs = cell(1, length(jsonmodel));
    linenum = cell(1, length(jsonmodel));
    sample = cell(1, length(jsonmodel));
    tagvalue = cell(1, length(jsonmodel));
    for i=1:length(jsonmodel)
        lhs{i} = jsonmodel{i}.lhs;
        rhs{i} = jsonmodel{i}.rhs;
        linenum{i} = jsonmodel{i}.line;
        if isfield(jsonmodel{i}, 'tags') && ...
                isfield(jsonmodel{i}.tags, 'name')
            tagvalue{i} = jsonmodel{i}.tags.('name');
        else
            tagvalue{i} = ['eq_line_no_' num2str(linenum{i})];
        end
        if isfield(jsonmodel{i}, 'tags')
            if isfield(jsonmodel{i}.tags, 'sample')
                tmp = strsplit(jsonmodel{i}.tags.sample, ':');
                sample{i} = dates(tmp{1}):dates(tmp{2});
            end
        else
            tagvalue{i} = ['eq_line_no_' num2str(linenum{i})];
        end
    end
    return
end

tagname = varargin{1};
tagvalue = varargin{2};

assert(ischar(tagname), 'Tag name must be a string');
assert(ischar(tagvalue) || iscell(tagvalue), 'Tag value must be a string or a cell string array');

if ischar(tagvalue)
    tagvalue = {tagvalue};
end

lhs = cell(1, length(tagvalue));
rhs = cell(1, length(tagvalue));
linenum = cell(1, length(tagvalue));
sample = cell(1, length(tagvalue));
idx2rm = [];
for j = 1:length(tagvalue)
    for i=1:length(jsonmodel)
        if isfield(jsonmodel{i}, 'tags') && ...
                isfield(jsonmodel{i}.tags, tagname) && ...
                strcmp(jsonmodel{i}.tags.(tagname), tagvalue{j})
            lhs{j} = jsonmodel{i}.lhs;
            rhs{j} = jsonmodel{i}.rhs;
            linenum{j} = jsonmodel{i}.line;
            if isfield(jsonmodel{i}.tags, 'sample')
                tmp = strsplit(jsonmodel{i}.tags.sample, ':');
                sample{j} = dates(tmp{1}):dates(tmp{2});
            end
            if ~any(cellfun(@isempty, lhs))
                return
            end
            break
        end
    end
    if isempty(rhs{j})
        warning(['getEquationsByTags: No equation tag found by the name of ''' tagvalue{j} ''''])
        idx2rm = [idx2rm j];
    end
end
if ~isempty(idx2rm)
    lhs(:,idx2rm) = [];
    rhs(:,idx2rm) = [];
    linenum(:,idx2rm) = [];
    sample(:,idx2rm) = [];
    tagvalue(:,idx2rm) = [];
end
end