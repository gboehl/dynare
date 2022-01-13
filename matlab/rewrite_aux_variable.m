function [variable, transformations] = rewrite_aux_variable(variable, M_)

% Copyright © 2021 Dynare Team
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

transformations = {};
ida = get_aux_variable_id(variable);
op = 0;
dl = 0;

while ida
    op = op+1;
    if isequal(M_.aux_vars(ida).type, 1)
        transformations(op,1) = {'lag'};
        transformations(op,2) = {M_.aux_vars(ida).orig_lead_lag-1};
        variable = M_.endo_names{M_.aux_vars(ida).orig_index};
    elseif isequal(M_.aux_vars(ida).type, 8)
        transformations(op, 1) = {'diff'};
        variable = M_.endo_names{M_.aux_vars(ida).orig_index};
    elseif isequal(M_.aux_vars(ida).type, 9)
        dl = dl+1;
        transformations(op, 1) = {'lag'};
        transformations(op, 2) = {dl};
        op = op-1;
        variable = M_.endo_names{M_.aux_vars(ida).orig_index};
    elseif isequal(M_.aux_vars(ida).type, 10)
        transformations(op, 1) = {M_.aux_vars(ida).unary_op};
        variable = M_.endo_names{M_.aux_vars(ida).orig_index};
    else
        error('This case is not implemented.')
    end
    ida = get_aux_variable_id(variable);
end