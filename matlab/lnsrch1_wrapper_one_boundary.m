function r = lnsrch1_wrapper_one_boundary(ya, ll_index, fname, y, x, params, steady_state, T, it_)
% wrapper for solve_one_boundary m-file when it is used with a dynamic
% model
%
% INPUTS
%   ya                  [vector]        The endogenous of the current block
%   ll_index            [vector]        M_.lead_lag_incidence(M_.maximum_endo_lag+1, :)
%   fname               [string]        name of the file containing the block
%                                       to simulate
%   y                   [vector]        Dynamic endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   steady_state        [vector]        steady state of the model
%   T                   [vector]        Temporary terms
% OUTPUTS
%   r                   [vector]        The residuals of the current block
%
% ALGORITHM
%   none.
%
% SPECIAL REQUIREMENTS
%   none.
%

% Copyright (C) 2009-2020 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licen

y2(nonzeros(ll_index)) = ya(find(ll_index));
[r, ~, g1]=feval(fname, y, x, params, steady_state, T, it_, false);
