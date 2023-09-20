function [steady_state, info, steady_state_changes_parameters] = check_steady_state_changes_parameters(M_,estim_params_,oo_,options_,steadystate_check_flag_vec)
% function [steady_state, info, steady_state_changes_parameters] = check_steady_state_changes_parameters(M_,estim_params_,oo_,options_,steadystate_check_flag_vec)
% -------------------------------------------------------------------------
% Check if steady-state solves static model and if it changes estimated parameters
% -------------------------------------------------------------------------
% INPUTS
%  o M_:                         [struct] information on the model
%  o estim_params_:              [struct] information on estimated parameters
%  o oo_:                        [struct] results
%  o options_:                   [struct] information on options
%  o steadystate_check_flag_vec: [vector] steadystate_check_flag for both checks (might be different in case of diffuse_filter)
% -------------------------------------------------------------------------
% OUTPUTS
%  o steady_state:                    [vector] steady state
%  o info:                            [scalar] 0 if steady state solves static model
%  o steady_state_changes_parameters: [logical] true if steady state file changes estimated parameters
% -------------------------------------------------------------------------
% This function is called by
%  o initial_estimation_checks.m
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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

value_parameter_change = 1.01; % value with which parameters are slightly changed.
steady_state_changes_parameters = false; % initialize

% check if steady state solves static model
[steady_state, new_steady_params, info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,steadystate_check_flag_vec(1));

% check whether steady state file changes estimated parameters
if isfield(estim_params_,'param_vals') && ~isempty(estim_params_.param_vals)
    old_steady_params = M_.params; % save initial parameters
    M_par_varied = M_; % store Model structure
    M_par_varied.params(estim_params_.param_vals(:,1)) = M_par_varied.params(estim_params_.param_vals(:,1))*value_parameter_change; % vary parameters
    [~, new_steady_params_2] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_par_varied,options_,steadystate_check_flag_vec(2));
    changed_par_indices = find((old_steady_params(estim_params_.param_vals(:,1))-new_steady_params(estim_params_.param_vals(:,1))) ...
                               | (M_par_varied.params(estim_params_.param_vals(:,1))-new_steady_params_2(estim_params_.param_vals(:,1))));
    if ~isempty(changed_par_indices)
        fprintf('\nThe steady state file internally changed the values of the following estimated parameters:\n')
        disp(char(M_.param_names(estim_params_.param_vals(changed_par_indices,1))))
        fprintf('This will override the parameter values and may lead to wrong results.\n')
        fprintf('Check whether this is really intended.\n')
        warning('The steady state file internally changes the values of the estimated parameters.')
        steady_state_changes_parameters = true;
    end
end