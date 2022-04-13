function b = isauxiliary(var, types)

% Returns true if var is an auxiliary variable.
%
% INPUTS
% - var       [string]    Name of the variable.
% - types     [integer]  vector of type of auxiliary variables.
%
% OUTPUTS
% - b         [logical]
%
% REMARKS
%
% Types for auxiliary variables are as follows:
%
%   -  0,    Lead on endogenous variable (substitute for endo leads >= 2)
%   -  1,    Lag on endogenous variable (ubstitute for endo lags >= 2)
%   -  2,    Lead on exogenous variable  (ubstitute for exo leads >= 1)
%   -  3,    Lag on exogenous variable (substitute for exo lags >= 1)
%   -  4     Expectation (substitute for Expectation Operator)
%   -  5,    Diff forward (substitute for the differentiate of a forward variable)
%   -  6,    Multipliers for FOC of Ramsey Problem
%   -  7,    currently not used, see SymbolTable.hh
%   -  8,    Variable for Diff operator
%   -  9,    Lag on Diff
%   - 10,    Variable created when diff was taken of unary operator (log, exp)
%   - 11,    Lead on Diff

% Copyright © 2018-2021 Dynare Team
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

b = false;

id = find(strcmp(var, M_.endo_names));

if isempty(id)
    return
end

if id<=M_.orig_endo_nbr
    return
else
    b = true;
    if nargin<2
        return
    end
end

auxinfo = M_.aux_vars(get_aux_variable_id(id));

if ~ismember(auxinfo.type, types)
    b = false;
end