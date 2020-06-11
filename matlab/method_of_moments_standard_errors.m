function [SEvalues, AVar] = method_of_moments_standard_errors(xparam, objective_function, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM, Wopt_flag)
% [SEvalues, AVar] = method_of_moments_standard_errors(xparam, objective_function, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM, Wopt_flag)
% -------------------------------------------------------------------------
% This function computes standard errors to the method of moments estimates
% Adapted from replication codes of
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
% =========================================================================
% INPUTS
%   o xparam:                   value of estimated parameters as returned by set_prior()
%   o objective_function        string of objective function, either method_of_moments_GMM.m or method_of_moments_SMM.m
%   o Bounds:                   structure containing parameter bounds
%   o DynareResults:            structure for results (oo_)
%   o EstimatedParameters:      structure describing the estimated_parameters (estim_params_)
%   o MatchedMoments:           structure containing information about selected moments to match in estimation (matched_moments_)
%   o Model                     structure describing the Model
%   o OptionsMoM:               structure information about all settings (specified by the user, preprocessor, and taken from global options_)
%   o Wopt_flag:                indicator whether the optimal weighting is actually used
% -------------------------------------------------------------------------
% OUTPUTS 
%   o SEvalues                  [nparam x 1] vector of standard errors
%   o AVar                      [nparam x nparam] asymptotic covariance matrix
% -------------------------------------------------------------------------
% This function is called by
%  o method_of_moments.m
% -------------------------------------------------------------------------
% This function calls:
%  o get_the_name
%  o get_error_message
%  o method_of_moments_GMM.m (objective function)
%  o method_of_moments_SMM.m (objective function)
%  o method_of_moments_optimal_weighting_matrix  
% =========================================================================
% Copyright (C) 2020 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Some dimensions
numMom      = size(DynareResults.mom.modelMoments,1);
dimParams   = size(xparam,1);
D           = zeros(numMom,dimParams);
epsValue    = OptionsMoM.dynatol.x;

for i=1:dimParams
    %Positive step
    xparam_eps_p      = xparam;
    xparam_eps_p(i,1) = xparam_eps_p(i) + epsValue;
    [~, info_p, exit_flag_p, DynareResults_p, ~, ~] = feval(objective_function, xparam_eps_p, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
    
    % Negative step
    xparam_eps_m      = xparam;
    xparam_eps_m(i,1) = xparam_eps_m(i) - epsValue;
    [~, info_m, exit_flag_m, DynareResults_m, ~, ~] = feval(objective_function, xparam_eps_m, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);

    % The Jacobian:
    if nnz(info_p)==0 && nnz(info_m)==0
        D(:,i) = (DynareResults_p.mom.modelMoments - DynareResults_m.mom.modelMoments)/(2*epsValue);
    else
        problpar = get_the_name(i,OptionsMoM.TeX, Model, EstimatedParameters, OptionsMoM);
        message_p = get_error_message(info_p, OptionsMoM);
        message_m = get_error_message(info_m, OptionsMoM);        
        
        warning('method_of_moments:info','Cannot compute the Jacobian for parameter %s - no standard errors available\n %s %s\nCheck your bounds and/or priors, or use a different optimizer.\n',problpar, message_p, message_m)
        AVar = NaN(length(xparam),length(xparam));
        SEvalues = NaN(length(xparam),1);
        return
    end
end

T = OptionsMoM.nobs; %Number of observations
if isfield(OptionsMoM,'variance_correction_factor')
    T = T*OptionsMoM.variance_correction_factor;
end

if Wopt_flag
    % We have the optimal weighting matrix    
    WW            = DynareResults.mom.Sw'*DynareResults.mom.Sw;
    AVar          = 1/T*((D'*WW*D)\eye(dimParams));
else
    % We do not have the optimal weighting matrix yet
    WWused        = DynareResults.mom.Sw'*DynareResults.mom.Sw;
    WWopt         = method_of_moments_optimal_weighting_matrix(DynareResults.mom.m_data, DynareResults.mom.modelMoments, OptionsMoM.bartlett_kernel_lag);
    S             = WWopt\eye(size(WWopt,1));
    AA            = (D'*WWused*D)\eye(dimParams);
    AVar          = 1/T*AA*D'*WWused*S*WWused*D*AA;
end

SEvalues   = sqrt(diag(AVar));