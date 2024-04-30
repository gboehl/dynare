function J_test = Jtest(xparam, objective_function, Q, model_moments, m_data, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% J_test = Jtest(xparam, objective_function, Q, model_moments, m_data, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% -------------------------------------------------------------------------
% Computes the J-test statistic and p-value given a GMM/SMM estimation
% -------------------------------------------------------------------------
% INPUTS
%  xparam:                [vector]           estimated parameter vector
%  objective_function:    [function handle]  objective function
%  Q:                     [scalar]           value of moments distance criterion
%  model_moments:         [vector]           model moments
%  m_data:                [matrix]           selected empirical moments at each point in time
%  data_moments:          [vector]           empirical moments
%  weighting_info:        [struct]           information on weighting matrix
%  options_mom_:          [struct]           options
%  M_:                    [struct]           model information
%  estim_params_:         [struct]           estimated parameters
%  bayestopt_:            [struct]           info on prior distributions
%  BoundsInfo:            [struct]           info bounds on parameters
%  dr:                    [struct]           reduced form model
%  endo_steady_state:     [vector]           steady state of endogenous variables (initval)
%  exo_steady_state:      [vector]           steady state of exogenous variables (initval)
%  exo_det_steady_state:  [vector]           steady state of deterministic exogenous variables (initval)
% -------------------------------------------------------------------------
% OUTPUT
%  J_test:              [struct]           results of J test
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%  o mom.objective_function
%  o mom.optimal_weighting_matrix
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


if options_mom_.mom.mom_nbr > length(xparam)
    % Get optimal weighting matrix for J test, if necessary
    if ~weighting_info.Woptflag
        W_opt = mom.optimal_weighting_matrix(m_data, model_moments, options_mom_.mom.bartlett_kernel_lag);
        weighting_info.Sw = chol(W_opt);
        fval = feval(objective_function, xparam, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    else
        fval = Q;
    end
    % Compute J statistic
    if strcmp(options_mom_.mom.mom_method,'SMM')
        variance_correction_factor = options_mom_.mom.variance_correction_factor;
    elseif strcmp(options_mom_.mom.mom_method,'GMM')
        variance_correction_factor = 1;
    end
    J_test.j_stat          = options_mom_.nobs*variance_correction_factor*fval/options_mom_.mom.weighting_matrix_scaling_factor;
    J_test.degrees_freedom = length(model_moments)-length(xparam);
    J_test.p_val           = 1-chi2cdf(J_test.j_stat, J_test.degrees_freedom);
    fprintf('\nValue of J-test statistic: %f\n',J_test.j_stat);
    fprintf('p-value of J-test statistic: %f\n',J_test.p_val);
else
    J_test=[];
end