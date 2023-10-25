function [endogenousvariables, success] = solve_stacked_linear_problem(endogenousvariables, exogenousvariables, steadystate_y, steadystate_x, M_, options_)

% Copyright © 2015-2023 Dynare Team
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

[options_, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, i_cols_0, i_cols_J0, dynamicmodel] = ...
    initialize_stacked_problem(endogenousvariables, options_, M_, steadystate_y);

ip = find(M_.lead_lag_incidence(1,:)');
ic = find(M_.lead_lag_incidence(2,:)');
in = find(M_.lead_lag_incidence(3,:)');

% Evaluate the Jacobian of the dynamic model at the deterministic steady state.
[d1,jacobian] = dynamicmodel(steadystate_y([ip; ic; in]), transpose(steadystate_x), M_.params, steadystate_y, 1);

% Check that the dynamic model was evaluated at the steady state.
if ~options_.steadystate.nocheck && max(abs(d1))>1e-12
    error('Jacobian is not evaluated at the steady state!')
end

nyp = nnz(M_.lead_lag_incidence(1,:));
ny0 = nnz(M_.lead_lag_incidence(2,:));
nyf = nnz(M_.lead_lag_incidence(3,:));
nd = nyp+ny0+nyf; % size of y (first argument passed to the dynamic file).
jexog = transpose(nd+(1:M_.exo_nbr));
jendo = transpose(1:nd);
z = bsxfun(@minus, z, steadystate_y);
x = bsxfun(@minus, exogenousvariables, steadystate_x');

[y, check, ~, ~, errorcode] = dynare_solve(@linear_perfect_foresight_problem, z(:), ...
                                           options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                           options_, ...
                                           jacobian, y0-steadystate_y, yT-steadystate_y, ...
                                           x, M_.params, steadystate_y, ...
                                           M_.maximum_lag, options_.periods, M_.endo_nbr, i_cols, ...
                                           i_cols_J1, i_cols_1, i_cols_T, i_cols_j, i_cols_0, i_cols_J0, ...
                                           jendo, jexog);

if all(imag(y)<.1*options_.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables = [y0 bsxfun(@plus,reshape(y,M_.endo_nbr,options_.periods), steadystate_y) yT];

success = ~check;

if ~success && options_.debug
    dprintf('solve_stacked_linear_problem: Nonlinear solver routine failed with errorcode=%i.', errorcode)
end
