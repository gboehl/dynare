function [pnames, enames, xnames, pid, eid, xid] = get_variables_and_parameters_in_equation(lhs, rhs, DynareModel)

% Returns the lists of parameters, endogenous variables and exogenous variables in an equation.
%
% INPUTS
% - lhs         [string]            Left hand side of an equation.
% - rhs         [string]            Right hand side of an equation.
% - DynareModel [struct]            Structure describing the current model (M_).
%
% OUTPUTS
% - pnames      [cell]              Cell of row char arrays (p elements), names of the parameters. 
% - enames      [cell]              Cell of row char arrays (n elements), names of the endogenous variables.
% - xnames      [cell]              Cell of row char arrays (m elements), names of the exogenous variables.
% - pid         [Integer]           p*1 vector of indices in M_.param_names for the listed parameters in params.
% - eid         [Integer]           n*1 vector of indices in M_.endo_names for the listed parameters in endogenous.
% - xid         [Integer]           m*1 vector of indices in M_.exo_names for the listed parameters in exogenous.

% Copyright (C) 2018-2019 Dynare Team
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

% Get the tokens in the rhs member of the equation.
rhs_ = strsplit(rhs,{'+','-','*','/','^', ...
                    'log(', 'log10(', 'ln(', 'exp(', ...
                    'sqrt(', 'abs(', 'sign(', ...
                    'sin(', 'cos(', 'tan(', 'asin(', 'acos(', 'atan(', ...
                    'min(', 'max(', ...
                    'normcdf(', 'normpdf(', 'erf(', ...
                    'diff(', 'adl(', '(', ')'});

% Filter out the numbers and punctuation.
rhs_(cellfun(@(x) all(isstrprop(x, 'digit')+isstrprop(x, 'punct')), rhs_)) = [];

% Get list of parameters.
pnames = DynareModel.param_names;
pnames = intersect(rhs_, pnames);

% Get list of endogenous variables.
enames = DynareModel.endo_names;
enames = intersect(rhs_, enames);

% Get list of exogenous variables
xnames = DynareModel.exo_names;
xnames = intersect(rhs_, xnames);

% Decide if we are dealing with a dynamic model. If so, the lhs variable
% already belongs to enames, we remove this variable from enames. 
id = find(strcmp(lhs, enames));
if ~isempty(id)
    enames(id) = [];
end

% Add lhs variable in first position of enames.
enames = [lhs; enames];

% Returns vector of indices for parameters endogenous and exogenous
% variables if required.
if nargout>3
    p = length(pnames);
    pid = zeros(p, 1);
    for i = 1:p
        pid(i) = find(strcmp(pnames{i}, DynareModel.param_names));
    end
    p = length(enames);
    eid = zeros(p, 1);
    for i = 1:p
        eid(i) = find(strcmp(enames{i}, DynareModel.endo_names));
    end
    p = length(xnames);
    xid = zeros(p, 1);
    for i = 1:p
        xid(i) = find(strcmp(xnames{i}, DynareModel.exo_names));
    end
end