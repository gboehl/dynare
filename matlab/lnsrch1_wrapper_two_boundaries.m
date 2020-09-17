function ra = lnsrch1_wrapper_two_boundaries(ya, fname, blk, y, y_index, x, ...
                                             params, steady_state, T, periods, ...
                                             y_size, M_)
% wrapper for solve_one_boundary m-file when it is used with a dynamic
% model
%
% INPUTS
%   ya                  [vector]        The endogenous of the current block
%   y_index             [vector of int] The index of the endogenous variables of
%                                       the block
%   fname               [string]        name of the dynamic file
%   blk                 [int]           block number
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   steady_state        [vector]        steady state of the model
%   T                   [matrix]        Temporary terms
%   periods             [int]           The number of periods
%   y_size              [int]           The number of endogenous variables
%                                       in the current block
%   M_                                  Model description structure
%
% OUTPUTS
%   ra                  [vector]        The residuals of the current block
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

%reshape the input arguments of the dynamic function
y(y_index, M_.maximum_lag+(1:periods)) = reshape(ya',length(y_index),periods);
ra = NaN(periods*y_size, 1);
for it_ = M_.maximum_lag+(1:periods)
    [ra((it_-M_.maximum_lag-1)*y_size+(1:y_size)), ~, ~, g1]=feval(fname, blk, dynvars_from_endo_simul(y, it_, M_), x, params, steady_state, T(:, it_), it_, false);
end
