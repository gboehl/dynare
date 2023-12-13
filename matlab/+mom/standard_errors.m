function [stderr_values, asympt_cov_mat] = standard_errors(xparam, objective_function, model_moments, model_moments_params_derivs, m_data, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [stderr_values, asympt_cov_mat] = standard_errors(xparam, objective_function, model_moments, model_moments_params_derivs, m_data, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% -------------------------------------------------------------------------
% This function computes standard errors to the method of moments estimates
% Adapted from replication codes of Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018):
% "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications",
% Review of Economic Studies, 85(1):1-49.
% -------------------------------------------------------------------------
% INPUTS
%  - xparam:               [vector]     value of estimated parameters as returned by set_prior()
%  - objective_function    [func]       function handle with string of objective function
%  - model_moments:        [vector]     model moments
%  - model_moments_params_derivs:  [matrix]  analytical jacobian of the model moments wrt estimated parameters (currently for GMM only)
%  - m_data                [matrix]     selected empirical moments at each point in time
%  - data_moments:         [vector]     data with moments/IRFs to match
%  - weighting_info:       [structure]  storing information on weighting matrices
%  - options_mom_:         [structure]  information about all settings (specified by the user, preprocessor, and taken from global options_)
%  - M_                    [structure]  model information
%  - estim_params_:        [structure]  information from estimated_params block
%  - bayestopt_:           [structure]  information on the prior distributions
%  - BoundsInfo:           [structure]  parameter bounds
%  - dr:                   [structure]  reduced form model
%  - endo_steady_state:    [vector]     steady state value for endogenous variables (initval)
%  - exo_steady_state:     [vector]     steady state value for exogenous variables (initval)
%  - exo_det_steady_state: [vector]     steady state value for exogenous deterministic variables (initval)
% -------------------------------------------------------------------------
% OUTPUTS 
%   o stderr_values        [nparam x 1] vector of standard errors
%   o asympt_cov_mat       [nparam x nparam] asymptotic covariance matrix
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run.m
% -------------------------------------------------------------------------
% This function calls:
%  o get_the_name
%  o get_error_message
%  o mom.objective_function
%  o mom.optimal_weighting_matrix
% -------------------------------------------------------------------------

% Copyright © 2020-2023 Dynare Team
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

% Some dimensions
num_mom      = size(model_moments,1);
dim_params   = size(xparam,1);
D            = zeros(num_mom,dim_params);
eps_value    = options_mom_.mom.se_tolx;

if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
    fprintf('\nComputing standard errors using analytical derivatives of moments\n');
    D = model_moments_params_derivs; % already computed in objective function via get_perturbation_params.m
    idx_nan = find(any(isnan(D)));
    if any(idx_nan)
        for i = idx_nan            
             fprintf('No standard errors available for parameter %s\n',get_the_name(i,options_mom_.TeX, M_, estim_params_, options_mom_.varobs))
        end        
        warning('There are NaN in the analytical Jacobian of Moments. Check your bounds and/or priors, or use a different optimizer.')
        asympt_cov_mat = NaN(length(xparam),length(xparam));
        stderr_values = NaN(length(xparam),1);
        return
    end
else    
    fprintf('\nComputing standard errors using numerical derivatives of moments\n');
    for i=1:dim_params
        % positive step
        xparam_eps_p      = xparam;
        xparam_eps_p(i,1) = xparam_eps_p(i) + eps_value;        
        [~, info_p, ~, ~, ~, ~, model_moments_p] = feval(objective_function, xparam_eps_p, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        % negative step
        xparam_eps_m      = xparam;
        xparam_eps_m(i,1) = xparam_eps_m(i) - eps_value;        
        [~, info_m, ~, ~, ~, ~, model_moments_m] = feval(objective_function, xparam_eps_m, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        % the Jacobian
        if nnz(info_p)==0 && nnz(info_m)==0
            D(:,i) = (model_moments_p - model_moments_m)/(2*eps_value);
        else
            problematic_parameter = get_the_name(i,options_mom_.TeX, M_, estim_params_, options_mom_.varobs);
            if info_p(1)==42
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s due to hitting the upper bound - no standard errors available.\n',problematic_parameter)
            else
                message_p = get_error_message(info_p, options_mom_);
            end
            if info_m(1)==41
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s due to hitting the lower bound - no standard errors available.\n',problematic_parameter)
            else
                message_m = get_error_message(info_m, options_mom_);
            end
            if info_m(1)~=41 && info_p(1)~=42
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s - no standard errors available\n %s %s\nCheck your priors or use a different optimizer.\n',problematic_parameter, message_p, message_m)
            end
            asympt_cov_mat = NaN(length(xparam),length(xparam));
            stderr_values = NaN(length(xparam),1);
            return
        end
    end
end
T = options_mom_.nobs;
if isfield(options_mom_,'variance_correction_factor')
    T = T*options_mom_.variance_correction_factor;
end
WW = weighting_info.Sw'*weighting_info.Sw;
if weighting_info.Woptflag
    % we already have the optimal weighting matrix
    asympt_cov_mat = 1/T*((D'*WW*D)\eye(dim_params));
else
    % we do not have the optimal weighting matrix yet
    WWopt      = mom.optimal_weighting_matrix(m_data, model_moments, options_mom_.mom.bartlett_kernel_lag);
    S          = WWopt\eye(size(WWopt,1));
    AA         = (D'*WW*D)\eye(dim_params);
    asympt_cov_mat = 1/T*AA*D'*WW*S*WW*D*AA;
end
stderr_values = sqrt(diag(asympt_cov_mat));