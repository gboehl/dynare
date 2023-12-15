function  [par, logpost, accepted, neval] = posterior_sampler_iteration(TargetFun,last_draw, last_posterior, sampler_options,varargin)

% function [par, logpost, accepted, neval] = posterior_sampler_iteration(TargetFun,last_draw, last_posterior, sampler_options,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_)
% posterior samplers
%
% INPUTS
%   TargetFun:              string storing the objective function (e.g. 'dsge_likelihood.m')
%   last_draw:              parameter vector in last iteration
%   last_posterior:         value of the posterior in last iteration
%   sampler_options:        posterior sampler options
%   dataset_:               the dataset after required transformation
%   dataset_info:           Various informations about the dataset (descriptive statistics and missing observations).
%   options_:               structure storing the options
%   M_:                     structure storing the model information
%   estim_params_:          structure storing information about estimated parameters
%   bayestopt_:             structure storing information about priors
%   mh_bounds:              structure containing prior bounds
%   oo_:                    structure storing the results
%
% OUTPUTS
%   par:                    last accepted parameter vector
%   logpost:                value of the posterior after current iteration
%   accepted:               share of proposed draws that were accepted
%   neval:                  number of evaluations (>1 only for slice)
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


posterior_sampling_method = sampler_options.posterior_sampling_method;
mh_bounds = sampler_options.bounds;

switch posterior_sampling_method
  case 'slice'

    [par, logpost, neval] = slice_sampler(TargetFun,last_draw, [mh_bounds.lb mh_bounds.ub], sampler_options,varargin{:});
    accepted = 1;
  case 'random_walk_metropolis_hastings'
    neval = 1;
    ProposalFun = sampler_options.proposal_distribution;
    proposal_covariance_Cholesky_decomposition = sampler_options.proposal_covariance_Cholesky_decomposition;
    n = sampler_options.n;

    par = feval(ProposalFun, last_draw, proposal_covariance_Cholesky_decomposition, n);
    if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
        try
            logpost = - feval(TargetFun, par(:),varargin{:});
        catch
            logpost = -inf;
        end
    else
        logpost = -inf;
    end
    r = logpost-last_posterior;
    if (logpost > -inf) && (log(rand) < r)
        accepted = 1;
    else
        accepted = 0;
        par = last_draw;
        logpost = last_posterior;
    end
  case 'tailored_random_block_metropolis_hastings'
    options_=varargin{3};
    bayestopt_=varargin{6};
    npar=length(last_draw);
    %% randomize indices for blocking in this iteration
    indices=randperm(npar)';
    blocks=[1; (1+cumsum((rand(length(indices)-1,1)>(1-sampler_options.new_block_probability))))];
    nblocks=blocks(end,1); %get number of blocks this iteration
    current_draw=last_draw'; %get starting point for current draw for updating
    blocked_draws_counter=0;
    accepted_draws_counter=0;
    for block_iter=1:nblocks
        blocked_draws_counter=blocked_draws_counter+1;
        nxopt=length(indices(blocks==block_iter,1)); %get size of current block
        par_start_current_block=current_draw(indices(blocks==block_iter,1));
        [xopt_current_block, fval, exitflag, hess_mat_optimizer, options_, Scale] = dynare_minimize_objective(@TaRB_optimizer_wrapper,par_start_current_block,sampler_options.mode_compute,options_,[mh_bounds.lb(indices(blocks==block_iter,1),1) mh_bounds.ub(indices(blocks==block_iter,1),1)],bayestopt_.name,bayestopt_,[],...
                                                          current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                                                          varargin{:}); %inputs for objective
        %% covariance for proposal density
        hessian_mat = reshape(hessian('TaRB_optimizer_wrapper',xopt_current_block, ...
                                      options_.gstep,...
                                      current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                                      varargin{:}),nxopt,nxopt);

        if any(any(isnan(hessian_mat))) || any(any(isinf(hessian_mat)))
            inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
        else
            inverse_hessian_mat=inv(hessian_mat); %get inverse Hessian
            if any(any((isnan(inverse_hessian_mat)))) || any(any((isinf(inverse_hessian_mat))))
                inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
            end
        end
        [proposal_covariance_Cholesky_decomposition_upper,negeigenvalues]=chol(inverse_hessian_mat);
        %if not positive definite, use generalized Cholesky of Eskow/Schnabel
        if negeigenvalues~=0
            proposal_covariance_Cholesky_decomposition_upper=chol_SE(inverse_hessian_mat,0);
        end
        proposal_covariance_Cholesky_decomposition_upper=proposal_covariance_Cholesky_decomposition_upper*diag(bayestopt_.jscale(indices(blocks==block_iter,1),:));
        %get proposal draw
        if strcmpi(sampler_options.proposal_distribution,'rand_multivariate_normal')
            n = nxopt;
        elseif strcmpi(sampler_options.proposal_distribution,'rand_multivariate_student')
            n = sampler_options.student_degrees_of_freedom;
        end

        proposed_par = feval(sampler_options.proposal_distribution, xopt_current_block', proposal_covariance_Cholesky_decomposition_upper, n);
        % check whether draw is valid and compute posterior
        if all( proposed_par(:) > mh_bounds.lb(indices(blocks==block_iter,1),:) ) && all( proposed_par(:) < mh_bounds.ub(indices(blocks==block_iter,1),:) )
            try
                logpost = - TaRB_optimizer_wrapper(proposed_par(:),...
                                  current_draw,indices(blocks==block_iter,1),TargetFun,...% inputs for wrapper
                                  varargin{:});
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end

        if (logpost > -inf)
            %get ratio of proposal densities, required because proposal depends
            %on current mode via Hessian and is thus not symmetric anymore
            if strcmpi(sampler_options.proposal_distribution,'rand_multivariate_normal')
                proposal_density_proposed_move_forward=multivariate_normal_pdf(proposed_par,xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
                proposal_density_proposed_move_backward=multivariate_normal_pdf(par_start_current_block',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
            elseif strcmpi(sampler_options.proposal_distribution,'rand_multivariate_student')
                proposal_density_proposed_move_forward=multivariate_student_pdf(proposed_par,xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
                proposal_density_proposed_move_backward=multivariate_student_pdf(par_start_current_block',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,n);
            end
            accprob=logpost-last_posterior+ log(proposal_density_proposed_move_backward)-log(proposal_density_proposed_move_forward); %Formula (6), Chib/Ramamurthy
            if  (log(rand) < accprob)
                current_draw(indices(blocks==block_iter,1))=proposed_par;
                last_posterior=logpost;
                accepted_draws_counter =accepted_draws_counter +1;
            else %no updating
                 %do nothing, keep old value
            end
        end
    end
    accepted=accepted_draws_counter/blocked_draws_counter;
    par = current_draw;
    neval=1;
    logpost = last_posterior; %make sure not a temporary draw is returned;
  case 'independent_metropolis_hastings'
    neval = 1;
    ProposalFun = sampler_options.proposal_distribution;
    ProposalDensity = sampler_options.ProposalDensity;
    proposal_covariance_Cholesky_decomposition = sampler_options.proposal_covariance_Cholesky_decomposition;
    n = sampler_options.n;
    xparam1 = sampler_options.xparam1';
    par = feval(ProposalFun, xparam1, proposal_covariance_Cholesky_decomposition, n);
    if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
        try
            logpost = - feval(TargetFun, par(:),varargin{:});
        catch
            logpost = -inf;
        end
    else
        logpost = -inf;
    end
    r = logpost - last_posterior + ...
        log(feval(ProposalDensity, last_draw, xparam1, proposal_covariance_Cholesky_decomposition, n)) - ...
        log(feval(ProposalDensity, par, xparam1, proposal_covariance_Cholesky_decomposition, n));
    if (logpost > -inf) && (log(rand) < r)
        accepted = 1;
    else
        accepted = 0;
        par = last_draw;
        logpost = last_posterior;
    end
end