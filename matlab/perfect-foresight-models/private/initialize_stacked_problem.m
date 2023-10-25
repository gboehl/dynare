function [options_, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, i_cols_0, i_cols_J0, dynamicmodel] = ...
    initialize_stacked_problem(endogenousvariables, options_, M_, steadystate_y)

% Sets up the stacked perfect foresight problem for use with dynare_solve.m
%
% INPUTS
% - endogenousvariables [double] N*T array, paths for the endogenous variables (initial guess).
% - options_            [struct] contains various options.
% - M_                  [struct] contains a description of the model.
% - steadystate_y       [double] N*1 array, steady state for the endogenous variables.
%
% OUTPUTS
% - options_            [struct] contains various options.
% - y0                  [double] N*1 array, initial conditions for the endogenous variables
% - yT                  [double] N*1 array, terminal conditions for the endogenous variables
% - z                   [double] T*M array, paths for the exogenous variables.
% - i_cols              [double] indices of variables appearing in M_.lead_lag_incidence
%                                and that need to be passed to _dynamic-file
% - i_cols_J1           [double] indices of contemporaneous and forward looking variables
%                                appearing in M_.lead_lag_incidence
% - i_cols_T            [double] columns of dynamic Jacobian related to
%                                contemporaneous and backward-looking
%                                variables (relevant in last period)
% - i_cols_j            [double] indices of variables in M_.lead_lag_incidence
%                                in dynamic Jacobian (relevant in intermediate periods)
% - i_cols_1            [double] indices of contemporaneous and forward looking variables in
%                                M_.lead_lag_incidence in dynamic Jacobian (relevant in first period)
% - i_cols_0            [double] indices of contemporaneous variables in M_.lead_lag_incidence in dynamic
%                                Jacobian (relevant in problems with periods=1)
% - i_cols_J0           [double] indices of contemporaneous variables appearing in M_.lead_lag_incidence (relevant in problems with periods=1)
% - dynamicmodel        [handle] function handle to _dynamic-file

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

periods = options_.periods;
if (options_.solve_algo == 10)
    if ~isfield(options_.lmmcp,'lb')
        [lb,ub,pfm.eq_index] = get_complementarity_conditions(M_,options_.ramsey_policy);
        if options_.linear_approximation
            lb = lb - steadystate_y;
            ub = ub - steadystate_y;
        end
        options_.lmmcp.lb = repmat(lb,periods,1);
        options_.lmmcp.ub = repmat(ub,periods,1);
    end
elseif (options_.solve_algo == 11)
    if ~isfield(options_.mcppath,'lb')
        [lb,ub,pfm.eq_index] = get_complementarity_conditions(M_,options_.ramsey_policy);
        if options_.linear_approximation
            lb = lb - steadystate_y;
            ub = ub - steadystate_y;
        end
        options_.mcppath.lb = repmat(lb,periods,1);
        options_.mcppath.ub = repmat(ub,periods,1);
    end
end

if M_.maximum_lag > 0
    y0 = endogenousvariables(:, M_.maximum_lag);
else
    y0 = NaN(M_.endo_nbr, 1);
end
if M_.maximum_lead > 0
    yT = endogenousvariables(:, M_.maximum_lag+periods+1);
else
    yT = NaN(M_.endo_nbr, 1);
end
z = endogenousvariables(:,M_.maximum_lag+(1:periods));
illi = M_.lead_lag_incidence';
[i_cols,~,i_cols_j] = find(illi(:));
if M_.maximum_endo_lag == 0
    i_cols = i_cols + M_.endo_nbr;
end
illi = illi(:,(1+M_.maximum_endo_lag):(1+M_.maximum_endo_lag+M_.maximum_endo_lead));
[i_cols_J1,~,i_cols_1] = find(illi(:));
i_cols_T = nonzeros(M_.lead_lag_incidence(1:(1+M_.maximum_endo_lag),:)');
if periods==1
    i_cols_0 = nonzeros(M_.lead_lag_incidence(1+M_.maximum_endo_lag,:)');
    i_cols_J0 = find(M_.lead_lag_incidence(1+M_.maximum_endo_lag,:)');
else
    i_cols_0 = [];
    i_cols_J0 = [];
end
dynamicmodel = str2func([M_.fname,'.dynamic']);
