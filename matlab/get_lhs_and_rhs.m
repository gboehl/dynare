function [lhs, rhs] = get_lhs_and_rhs(eqname, DynareModel, original)

% Returns the left and right handsides of an equation.
%
% INPUTS
% - lhs         [string]            Left hand side of the equation.
% - rhs         [string]            Right hand side of the equation.
% - DynareModel [struct]            Structure describing the current model (M_).
%
% OUTPUTS
% - eqname      [string]            Name of the equation.
%
% SPECIAL REQUIREMENTS
%  The user must have attached names to the equations using equation
%  tags. Each equation in the model block must be preceeded with a
%  tag (see the reference manual). For instance, we should have
%  something as:
%
%      [name='Phillips curve']
%      pi = beta*pi(1) + slope*y + lam;

% Copyright (C) 2018 Dynare Team
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

if nargin<3
    original = false;
end

% Get the equation from the JSON output.
if original
    jsonfil = loadjson([DynareModel.fname '_original.json']);
else
    jsonfil = loadjson([DynareModel.fname '.json']);
end
jsonmod = jsonfil.model;
jsoneqn = getEquationsByTags(jsonmod, 'name', eqname);

% Get the lhs and rhs members of the selected equation.
lhs = jsoneqn{1}.lhs;
rhs = jsoneqn{1}.rhs;