function [forecast, error_flag] = forecast(options_,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,forecast_horizon)
% [forecast, error_flag] = forecast(options_,M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,forecast_horizon)
% Occbin forecasts
%
% INPUTS
% - options_                [structure]     Matlab's structure describing the current options
% - M_                      [structure]     Matlab's structure describing the model
% - dr_in                   [structure]     model information structure
% - endo_steady_state       [double]        steady state value for endogenous variables
% - exo_steady_state        [double]        steady state value for exogenous variables
% - exo_det_steady_state    [double]        steady state value for exogenous deterministic variables                                    
% - forecast_horizon        [integer]       forecast horizon
%
% OUTPUTS
% - forecast                [structure]     forecast results
% - error_flag              [integer]       error code
%
% SPECIAL REQUIREMENTS
%   none.

% Copyright © 2022-2023 Dynare Team
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

opts = options_.occbin.forecast;

options_.occbin.simul.maxit = opts.maxit;
options_.occbin.simul.check_ahead_periods = opts.check_ahead_periods;
options_.occbin.simul.periods = forecast_horizon;
shocks_input = opts.SHOCKS0;

if ~isempty(shocks_input)
    n_shocks=size(shocks_input,1);
    if iscell(shocks_input)
        inds=NaN(n_shocks,1);
        periods=length(shocks_input{1}{2});
        shock_mat=NaN(n_shocks,periods);
        for j=1:n_shocks
            exo_pos=strmatch(shocks_input{j}{1},M_.exo_names,'exact');
            if isempty(exo_pos)
        	    error('occbin.forecast: unknown exogenous shock %s',shocks_input{j}{1})
            else
                inds(j)=exo_pos;
            end
            if length(shocks_input{j}{2})~=periods
        	    error('occbin.forecast: unknown exogenous shock %s',shocks_input{j}{1})
            else
                shock_mat(j,:)=shocks_input{j}{2};
            end
        end
    elseif isreal(shocks_input)
        shock_mat=shocks_input;
        inds = (1:M_.exo_nbr)';
    end
end

options_.occbin.simul.endo_init = M_.endo_histval(:,1)-endo_steady_state; %initial condition
options_.occbin.simul.init_regime = opts.frcst_regimes;
options_.occbin.simul.init_binding_indicator = [];

shocks_base = zeros(forecast_horizon,M_.exo_nbr);
if ~isempty(shocks_input)
    for j=1:n_shocks
        shocks_base(:,inds(j))=shock_mat(j,:);
    end
end

if opts.replic
    h = dyn_waitbar(0,'Please wait occbin forecast replic ...');
    ishock = find(sqrt(diag((M_.Sigma_e))));
    options_.occbin.simul.exo_pos=ishock;
    effective_exo_nbr=  length(ishock);
    effective_Sigma_e = M_.Sigma_e(ishock,ishock);  % does not take heteroskedastic shocks into account
    [U,S] = svd(effective_Sigma_e);
    % draw random shocks
    if opts.qmc
        opts.replic =2^(round(log2(opts.replic+1)))-1;
        SHOCKS_add = qmc_sequence(forecast_horizon*effective_exo_nbr, int64(1), 1, opts.replic);
    else
        SHOCKS_add = randn(forecast_horizon*effective_exo_nbr,opts.replic);
    end
    SHOCKS_add=reshape(SHOCKS_add,effective_exo_nbr,forecast_horizon,opts.replic);
    z.linear=NaN(forecast_horizon,M_.endo_nbr,opts.replic);
    z.piecewise=NaN(forecast_horizon,M_.endo_nbr,opts.replic);
    error_flag=true(opts.replic,1);
    simul_SHOCKS=NaN(forecast_horizon,M_.exo_nbr,opts.replic);
    for iter=1:opts.replic
        options_.occbin.simul.SHOCKS = shocks_base+transpose(U*sqrt(S)*SHOCKS_add(:,:,iter));
        options_.occbin.simul.waitbar=0;
        [~, out] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
        error_flag(iter)=out.error_flag;
        if ~error_flag(iter)
            z.linear(:,:,iter)=out.linear;
            z.piecewise(:,:,iter)=out.piecewise;
            frcst_regime_history(iter,:)=out.regime_history;
            error_flag(iter)=out.error_flag;
            simul_SHOCKS(:,:,iter) = shocks_base;
        else
            if options_.debug
                save('Occbin_forecast_debug','simul_SHOCKS','z','iter','frcst_regime_history','error_flag','out','shocks_base')
            end
        end
        dyn_waitbar(iter/opts.replic,h,['OccBin MC forecast replic ',int2str(iter),'/',int2str(opts.replic)])
    end
    dyn_waitbar_close(h);
    if options_.debug
         save('Occbin_forecast_debug','simul_SHOCKS','z','iter','frcst_regime_history','error_flag')
    end
    inx=find(error_flag==0);
    z.linear=z.linear(:,:,inx);
    z.piecewise=z.piecewise(:,:,inx);
    z.min.piecewise = min(z.piecewise,[],3);
    z.max.piecewise = max(z.piecewise,[],3);
    z.min.linear = min(z.linear,[],3);
    z.max.linear = max(z.linear,[],3);
    
    field_names={'linear','piecewise'};
    post_mean=NaN(forecast_horizon,1);
    post_median=NaN(forecast_horizon,1);
    post_var=NaN(forecast_horizon,1);
    hpd_interval=NaN(forecast_horizon,2);
    post_deciles=NaN(forecast_horizon,9);
    for field_iter=1:2
        for i=1:M_.endo_nbr
            for j=1:forecast_horizon
                [post_mean(j,1), post_median(j,1), post_var(j,1), hpd_interval(j,:), post_deciles(j,:)] = posterior_moments(squeeze(z.(field_names{field_iter})(j,i,:)),options_.forecasts.conf_sig);
            end
            forecast.(field_names{field_iter}).Mean.(M_.endo_names{i})=post_mean;
            forecast.(field_names{field_iter}).Median.(M_.endo_names{i})=post_median;
            forecast.(field_names{field_iter}).Var.(M_.endo_names{i})=post_var;
            forecast.(field_names{field_iter}).HPDinf.(M_.endo_names{i})=hpd_interval(:,1);
            forecast.(field_names{field_iter}).HPDsup.(M_.endo_names{i})=hpd_interval(:,2);
            forecast.(field_names{field_iter}).Deciles.(M_.endo_names{i})=post_deciles;
            forecast.(field_names{field_iter}).Min.(M_.endo_names{i})=z.min.(field_names{field_iter})(:,i);
            forecast.(field_names{field_iter}).Max.(M_.endo_names{i})=z.max.(field_names{field_iter})(:,i);
        end
    end
else
    options_.occbin.simul.irfshock = M_.exo_names;
    options_.occbin.simul.SHOCKS = shocks_base;
    [~, out] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
    error_flag=out.error_flag;
    if out.error_flag
        fprintf('occbin.forecast: forecast simulation failed.')
        return;
    end

    frcst_regime_history=out.regime_history;
    error_flag=out.error_flag;
    for i=1:M_.endo_nbr
        forecast.linear.Mean.(M_.endo_names{i})= out.linear(:,i);
        forecast.piecewise.Mean.(M_.endo_names{i})= out.piecewise(:,i);
    end
end

forecast.regimes=frcst_regime_history;