function measure = measurement_equations(StateVectors,ReducedForm,ThreadsOptions, options_, M_)

% Copyright © 2013-2022 Dynare Team
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

order = options_.order;
mf1 = ReducedForm.mf1;
if ReducedForm.use_k_order_solver
    dr = ReducedForm.dr;
    udr = ReducedForm.udr;
else
    ghx  = ReducedForm.ghx(mf1,:);
    ghu  = ReducedForm.ghu(mf1,:);
    ghxx = ReducedForm.ghxx(mf1,:);
    ghuu = ReducedForm.ghuu(mf1,:);
    ghxu = ReducedForm.ghxu(mf1,:);
    ghs2 = ReducedForm.ghs2(mf1,:);
    if order == 3
        ghxxx = ReducedForm.ghxxx(mf1,:);
        ghuuu = ReducedForm.ghuuu(mf1,:);
        ghxxu = ReducedForm.ghxxu(mf1,:);
        ghxuu = ReducedForm.ghxuu(mf1,:);
        ghxss = ReducedForm.ghxss(mf1,:);
        ghuss = ReducedForm.ghuss(mf1,:);
    end
end
steadystate = ReducedForm.steadystate(mf1,:);
constant = ReducedForm.constant(mf1,:);
state_variables_steady_state = ReducedForm.state_variables_steady_state;
number_of_structural_innovations = length(ReducedForm.Q);
yhat = bsxfun(@minus, StateVectors, state_variables_steady_state);
if ReducedForm.use_k_order_solver
    tmp = local_state_space_iteration_k(yhat, zeros(number_of_structural_innovations, size(yhat,2)), dr, M_, options_, udr);
    measure = tmp(mf1,:);
else
    if order == 2
        measure = local_state_space_iteration_2(yhat, zeros(number_of_structural_innovations, size(yhat,2)), ghx, ghu, constant, ghxx, ghuu, ghxu, ThreadsOptions.local_state_space_iteration_2);
    elseif order == 3
        measure = local_state_space_iteration_3(yhat, zeros(number_of_structural_innovations, size(yhat,2)), ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, false);
    else
        error('Order > 3: use_k_order_solver should be set to true');
    end
end
