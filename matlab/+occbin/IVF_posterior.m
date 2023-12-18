function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,M_,options_,bayestopt_,dr, atT, innov] = IVF_posterior(xparam1,...
    dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,M_,options_,bayestopt_,dr, atT, innov] = IVF_posterior(xparam1,...
%     dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Computes Likelihood with inversion filter
%
% INPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - dataset_            [structure]     dataset after transformations
% - dataset_info        [structure]     storing informations about the
%                                       sample; not used but required for interface
% - options_            [structure]     Matlab's structure describing the current options
% - M_                  [structure]     Matlab's structure describing the model
% - estim_params_       [structure]     characterizing parameters to be estimated
% - bayestopt_          [structure]     describing the priors
% - BoundsInfo          [structure]     containing prior bounds
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state [vector]       steady state value for exogenous deterministic variables
%
% OUTPUTS
% - fval                    [double]        scalar, value of the likelihood or posterior kernel.
% - info                    [integer]       4×1 vector, informations resolution of the model and evaluation of the likelihood.
% - exit_flag               [integer]       scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - DLIK                    [double]        Empty array.
% - Hess                    [double]        Empty array.
% - SteadyState             [double]        Empty array.
% - trend                   [double]        Empty array.
% - M_                      [struct]        Updated M_ structure described in INPUTS section.
% - options_                [struct]        Updated options_ structure described in INPUTS section.
% - bayestopt_              [struct]        See INPUTS section.
% - dr                      [structure]     Reduced form model.
% - atT                     [double]        (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
% - innov                   [double]        (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).

% Copyright © 2021-2023 Dynare Team
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


DLIK=[];
Hess=[];
trend_coeff = [];
obs = dataset_.data;
obs_list = options_.varobs(:);
exit_flag   = 1;


if size(xparam1,1)<size(xparam1,2)
    xparam1=xparam1';
end


%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

if ~isempty(xparam1)
    M_ = set_all_parameters(xparam1,estim_params_,M_);
    [fval,info,exit_flag]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, BoundsInfo);
    if info(1)
        return
    end
end

err_index=options_.occbin.likelihood.IVF_shock_observable_mapping; % err_index= find(diag(M_.Sigma_e)~=0);
COVMAT1 = M_.Sigma_e(err_index,err_index);

% Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
[~,~,SteadyState,info,dr, M_.params] = dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end
end

%---------------------------------------------
% Calculate likelihood
%---------------------------------------------
sample_length = size(obs,1);
filtered_errs_init = zeros(sample_length,sum(err_index));

[filtered_errs, resids, Emat, stateval, info] = occbin.IVF_core(M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,err_index,filtered_errs_init,obs_list,obs);
if info(1)
    fval = Inf;
    exit_flag = 0;
    atT=NaN(size(stateval(:,dr.order_var)'));
    innov=NaN(M_.exo_nbr,sample_length);
    return
else
    atT = stateval(:,dr.order_var)';
    innov = zeros(M_.exo_nbr,sample_length);
    innov(diag(M_.Sigma_e)~=0,:)=filtered_errs';
end
nobs=size(filtered_errs,1);


%-------------------------------------
% Calculate the selection matrix
%-------------------------------------

iobs=options_.varobs_id';

likei=NaN(nobs,1);
if ~any(any(isnan(obs)))
    log_det_COV=log(det(COVMAT1));
    for t = 1:nobs
        Gmat1  = Emat(iobs,:,t);
        log_det_jacobian = log_det_COV + 2*log(abs(det(Gmat1)));
        trace_term = filtered_errs(t,:)/(COVMAT1)*filtered_errs(t,:)';
        
        likei(t,1) = log_det_jacobian + trace_term;
    end
    likei=likei+length(iobs)*log(2*pi);
else
    for t = 1:nobs
        inan=~isnan(obs(t,:));
        Gmat1  = Emat(iobs(inan),inan,t);
        log_det_jacobian = log(det(COVMAT1(inan,inan))) + 2*log(abs(det(Gmat1)));
        trace_term = filtered_errs(t,inan)/(COVMAT1(inan,inan))*filtered_errs(t,inan)';
        
        likei(t,1) = log_det_jacobian + trace_term + length(iobs(inan))*log(2*pi);
    end    
end

like = 0.5*sum(likei(options_.presample+1:end));

if isinf(like) 
    fval = Inf;
    info(1) = 301;
    info(4) = 1000;
    exit_flag = 0;
    return
elseif isnan(like)
    fval = Inf;
    info(1) = 302;
    info(4) = 1000;
    exit_flag = 0;
    return
end

maxresid = max(abs(resids(:)));
if maxresid>1e-3
    disp_verbose('Penalize failure of residuals to be zero',options_.verbosity)    
    fval = Inf;
    info(1) = 303;
    info(4) = sum(resids(:).^2);
    exit_flag = 0;
    return
end

if ~isempty(xparam1)
    prior = -priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
else
    prior = 0;
end
if prior == Inf
    % If parameters outside prior bound, minus prior density is very large
    fval = Inf;
    info(4) = 1000;
    exit_flag = 0;
    return
end

%---------------------------------------------
% Calculate posterior
%---------------------------------------------

% remember that the likelihood has already been multiplied by -1
% hence, posterior is -1 times the log of the prior
fval = like+prior;
