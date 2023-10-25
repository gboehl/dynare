function [lhs, rhs, json] = get_lhs_and_rhs(eqname, M_, original, json)
% [lhs, rhs, json] = get_lhs_and_rhs(eqname, M_, original, json)
% Returns the left and right handsides of an equation.
%
% INPUTS
% - eqname       [char]            Name of the equation.
% - M_           [struct]          Structure describing the current model.
% - original     [logical]         fetch equation in modfile-original.json or modfile.json
% - json         [char]            content of the JSON file
%
% OUTPUTS
% - lhs          [char]            Left hand side of the equation.
% - rhs          [char]            Right hand side of the equation.

%
%
% SPECIAL REQUIREMENTS
%  The user must have attached names to the equations using equation
%  tags. Each equation in the model block must be preceeded with a
%  tag (see the reference manual). For instance, we should have
%  something as:
%
%      [name='Phillips curve']
%      pi = beta*pi(1) + slope*y + lam;

% Copyright © 2018-2023 Dynare Team
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

if nargin<3 || isempty(original)
    original = false;
end

% Load JSON file if nargin<4
if nargin<4
    if original
        json = loadjson_([M_.fname filesep() 'model' filesep() 'json' filesep() 'modfile-original.json']);
    else
        json = loadjson_([M_.fname filesep() 'model' filesep() 'json' filesep() 'modfile.json']);
    end
end

% Load model.
jsonmod = json.model;
if isstruct(jsonmod)
    jsonmod = {jsonmod};
end

% Load equation.
jsoneqn = getEquationsByTags(jsonmod, 'name', eqname);

% Get the lhs and rhs members of the selected equation.
lhs = jsoneqn{1}.lhs;
rhs = jsoneqn{1}.rhs;