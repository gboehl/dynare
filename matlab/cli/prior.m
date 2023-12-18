function varargout = prior(varargin)
% varargout = prior(varargin)
% Computes various prior statistics and display them in the command window.
%
% INPUTS
%   'table', 'moments', 'optimize', 'simulate', 'plot', 'moments(distribution)'
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

if isempty(varargin) || ( isequal(length(varargin), 1) && isequal(varargin{1},'help'))
    skipline()
    disp('Possible options are:')
    disp(' + table                 Prints a table describing the priors.')
    disp(' + moments               Computes and displays moments of the endogenous variables at the prior mode.')
    disp(' + optimize              Optimizes the prior density (starting from a random initial guess).')
    disp(' + simulate              Computes the effective prior mass (using a Monte-Carlo).')
    disp(' + plot                  Plots the marginal prior densities.')
    disp(' + moments(distribution) Print tables describing the implied prior for the first and second order unconditional')
    disp('                         moments of all the endogenous variables.')
    skipline()
    return
end

global options_ M_ estim_params_ bayestopt_ oo_

donesomething = false;

if ~isbayes(estim_params_)
    warning('No prior detected!')
    return
end

if (size(estim_params_.var_endo,1) || size(estim_params_.corrn,1))
    % Prior over measurement errors are defined...
    if ((isfield(options_,'varobs') && isempty(options_.varobs)) || ~isfield(options_,'varobs'))
        % ... But the list of observed variabled is not yet defined.
        warning('Prior detected on measurement erros, but no list of observed variables (varobs is missing)!')
        return
    end
end

% Fill or update bayestopt_ structure
[~, estim_params_, ~, lb, ub, M_local] = set_prior(estim_params_, M_, options_);
% Set restricted state space
options_plot_priors_old=options_.plot_priors;
options_.plot_priors=0;
[~,~,~,~, M_, options_, oo_, estim_params_, BayesOptions] = ...
    dynare_estimation_init(M_.endo_names, M_.fname, 1, M_, options_, oo_, estim_params_, bayestopt_);
options_.plot_priors=options_plot_priors_old;


% Temporarly change qz_criterium option value
changed_qz_criterium_flag  = 0;
if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-9;
    changed_qz_criterium_flag  = 1;
end

M_local.dname = M_local.fname;

% Temporarly set options_.order equal to one
order = options_.order;
options_.order = 1;

if ismember('plot', varargin)
    plot_priors(BayesOptions, M_local, estim_params_, options_)
    donesomething = true;
end

if ismember('table', varargin)
    print_table_prior(lb, ub, options_, M_local, BayesOptions, estim_params_);
    donesomething = true;
end

if ismember('simulate', varargin) % Prior simulations (BK).
    if ismember('moments(distribution)', varargin)
        results = prior_sampler(1, M_local, BayesOptions, options_, oo_, estim_params_);
    else
        results = prior_sampler(0, M_local, BayesOptions, options_, oo_, estim_params_);
    end
    % Display prior mass info
    skipline(2)
    disp(['Prior mass = ' num2str(results.prior.mass)])
    disp(['BK indeterminacy share                = ' num2str(results.bk.indeterminacy_share)])
    disp(['BK unstability share                  = ' num2str(results.bk.unstability_share)])
    disp(['BK singularity share                  = ' num2str(results.bk.singularity_share)])
    disp(['Complex jacobian share                = ' num2str(results.jacobian.problem_share)])
    disp(['mjdgges crash share                   = ' num2str(results.dll.problem_share)])
    disp(['Steady state problem share            = ' num2str(results.ss.problem_share)])
    disp(['Complex steady state share            = ' num2str(results.ss.complex_share)])
    disp(['Endogenous prior violation share      = ' num2str(results.endogenous_prior_violation_share)])
    if options_.loglinear
        disp(['Nonpositive steady state share        = ' num2str(results.ss.nonpositive_share)])
    end
    disp(['Analytical steady state problem share = ' num2str(results.ass.problem_share)])
    skipline(2)
    donesomething = true;
end

if ismember('optimize', varargin) % Prior optimization.
    optimize_prior(options_, M_local, oo_, BayesOptions, estim_params_);
    donesomething = true;
end

if ismember('moments', varargin) % Prior simulations (2nd order moments).
                                 % Set estimated parameters to the prior mode...
    xparam1 = BayesOptions.p5;
    % ... Except for uniform priors (use the prior mean)!
    k = find(isnan(xparam1));
    xparam1(k) = BayesOptions.p1(k);
    % Update vector of parameters and covariance matrices
    M_local = set_all_parameters(xparam1, estim_params_, M_local);
    % Check model.
    check_model(M_local);
    % Compute state space representation of the model.
    oo__ = oo_;
    oo__.dr = set_state_space(oo__.dr, M_local);
    % Solve model
    [T,R,~,info,oo__.dr, M_local.params] = dynare_resolve(M_local , options_ , oo__.dr, oo__.steady_state, oo__.exo_steady_state, oo__.exo_det_steady_state,'restrict');
    if ~info(1)
        info=endogenous_prior_restrictions(T,R,M_local , options__ , oo__.dr,oo__.steady_state,oo__.exo_steady_state,oo__.exo_det_steady_state);
    end
    if info
        skipline()
        message = get_error_message(info,options_);
        fprintf('Cannot solve the model on the prior mode (info = %d, %s)\n', info(1), message);
        skipline()
        return
    end
    % Compute and display second order moments
    oo__ = disp_th_moments(oo__.dr, [], M_local, options__, oo__);
    skipline(2)
    donesomething = true;
end

if ismember('moments(distribution)', varargin) % Prior simulations (BK).
    if ~ismember('simulate', varargin)
        results = prior_sampler(1, M_local, BayesOptions, options_, oo_, estim_params_);
    end
    priorpath = [M_local.dname filesep() 'prior' filesep() 'draws' filesep()];
    list_of_files = dir([priorpath 'prior_draws*']);
    FirstOrderMoments = NaN(M_local.orig_endo_nbr, options_.prior_mc);
    SecondOrderMoments = NaN(M_local.orig_endo_nbr, M_local.orig_endo_nbr, options_.prior_mc);
    iter = 1;
    noprint = options_.noprint;
    options_.noprint = 1;
    for i=1:length(list_of_files)
        tmp = load([priorpath list_of_files(i).name]);
        for j = 1:size(tmp.pdraws, 1)
            if ~tmp.pdraws{j,2}
                dr = tmp.pdraws{j,3};
                oo__ = oo_;
                oo__.dr = dr;
                M_local=set_parameters_locally(M_local,tmp.pdraws{j,1});% Needed to update the covariance matrix of the state innovations.
                oo__ = disp_th_moments(oo__.dr, [], M_local, options_, oo__);
                FirstOrderMoments(:,iter) = oo__.mean;
                SecondOrderMoments(:,:,iter) = oo__.var;
                iter = iter+1;
            end
        end
    end
    save([M_.dname filesep() 'prior' filesep() M_.fname '_endogenous_variables_prior_draws.mat'], 'FirstOrderMoments', 'SecondOrderMoments')
    skipline(2)
    options_.noprint = noprint;
    % First order moments
    FirstOrderMoments = FirstOrderMoments(:,1:iter-1);
    SecondOrderMoments = SecondOrderMoments(:,:,1:iter-1);
    PriorExpectationOfFirstOrderMoments = mean(FirstOrderMoments, 2);
    PriorVarianceOfFirstOrderMoments = ...
        mean(bsxfun(@minus, FirstOrderMoments, PriorExpectationOfFirstOrderMoments).^2, 2);
    % Second order moments
    PriorExpectationOfSecondOrderMoments = mean(SecondOrderMoments, 3);
    PriorVarianceOfSecondOrderMoments = ...
        mean(bsxfun(@minus, SecondOrderMoments, PriorExpectationOfSecondOrderMoments).^2, 3);
    % Display first and second order moments implied priors (expectation and variance)
    print_moments_implied_prior(M_, PriorExpectationOfFirstOrderMoments, ...
                                PriorVarianceOfFirstOrderMoments, ...
                                PriorExpectationOfSecondOrderMoments, ...
                                PriorVarianceOfSecondOrderMoments);
    donesomething = true;
end

if changed_qz_criterium_flag
    options_.qz_criterium = [];
end

options_.order = order;

if ~donesomething
    error('prior: Unexpected arguments!')
end