function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,Model,DynareOptions,BayesInfo,DynareResults, atT, innov] = IVF_posterior(xparam1,...
    dataset_,obs_info,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults)
% function [fval,info,exit_flag,DLIK,Hess,SteadyState,trend_coeff,Model,DynareOptions,BayesInfo,DynareResults, atT, innov] = IVF_posterior(xparam1,...
%     dataset_,obs_info,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults)
% Computes Likelihood with inversion filter
%
% INPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - dataset_            [structure]     dataset after transformations
% - DynareOptions       [structure]     Matlab's structure describing the current options (options_).
% - Model               [structure]     Matlab's structure describing the model (M_).
% - EstimatedParameters [structure]     characterizing parameters to be estimated
% - BayesInfo           [structure]     describing the priors
% - BoundsInfo          [structure]     containing prior bounds
% - DynareResults       [structure]     Matlab's structure containing the results (oo_).
%
% OUTPUTS
% - fval                    [double]        scalar, value of the likelihood or posterior kernel.
% - info                    [integer]       4×1 vector, informations resolution of the model and evaluation of the likelihood.
% - exit_flag               [integer]       scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - DLIK                    [double]        Empty array.
% - Hess                    [double]        Empty array.
% - SteadyState             [double]        Empty array.
% - trend                   [double]        Empty array.
% - Model                   [struct]        Updated Model structure described in INPUTS section.
% - DynareOptions           [struct]        Updated DynareOptions structure described in INPUTS section.
% - BayesInfo               [struct]        See INPUTS section.
% - DynareResults           [struct]        Updated DynareResults structure described in INPUTS section.
% - atT                     [double]        (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
% - innov                   [double]        (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).

% Copyright (C) 2021 Dynare Team
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
obs_list = DynareOptions.varobs(:);
exit_flag   = 1;


if size(xparam1,1)<size(xparam1,2)
    xparam1=xparam1';
end


%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

if ~isempty(xparam1)
    Model = set_all_parameters(xparam1,EstimatedParameters,Model);
    [fval,info,exit_flag,Q,H]=check_bounds_and_definiteness_estimation(xparam1, Model, EstimatedParameters, BoundsInfo);
    if info(1)
        return
    end
end

err_index=DynareOptions.occbin.likelihood.IVF_shock_observable_mapping; % err_index= find(diag(Model.Sigma_e)~=0);
COVMAT1 = Model.Sigma_e(err_index,err_index);

% Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
[T,R,SteadyState,info,Model,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

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

[filtered_errs, resids, Emat, stateval, info] = occbin.IVF_core(Model,DynareResults,DynareOptions,err_index,filtered_errs_init,obs_list,obs);
if info(1)
    fval = Inf;
    exit_flag = 0;
    return
end
nobs=size(filtered_errs,1);


%-------------------------------------
% Calculate the selection matrix
%-------------------------------------

iobs=DynareOptions.varobs_id';

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

like = 0.5*sum(likei(DynareOptions.presample+1:end));

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
    prior = -priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
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
atT = stateval(:,DynareResults.dr.order_var)';
innov = zeros(Model.exo_nbr,sample_length);
innov(diag(Model.Sigma_e)~=0,:)=filtered_errs';
updated_variables = atT*nan;
BayesInfo.mf = BayesInfo.smoother_var_list(BayesInfo.smoother_mf);


initDynareOptions=DynareOptions;
DynareOptions.nk=[]; %unset options_.nk and reset it later
[DynareResults]=store_smoother_results(Model,DynareResults,DynareOptions,BayesInfo,dataset_,obs_info,atT,innov,[],updated_variables,DynareResults.dr.ys,zeros(length(DynareOptions.varobs_id),1));
DynareOptions=initDynareOptions;
