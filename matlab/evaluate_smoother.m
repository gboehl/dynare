function [oo_,M_,options_,bayestopt_,Smoothed_variables_declaration_order_deviation_form,initial_date]=evaluate_smoother(parameters,var_list,M_,oo_,options_,bayestopt_,estim_params_)
% Evaluate the smoother at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean','mle_mode','calibration') or a vector of values for
%                  the (estimated) parameters of the model.
%    o var_list    subset of endogenous variables
%    o M_          [structure]  Definition of the model
%    o oo_         [structure]  Storage of results
%    o options_    [structure]  Options
%    o bayestopt_  [structure]  describing the priors
%    o estim_params_ [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    o oo       [structure]  results:
%                              - SmoothedVariables
%                              - SmoothedShocks
%                              - FilteredVariablesShockDecomposition
%                              - UpdatedVariables
%                              - FilteredVariables
%                              - SmoothedMeasurementErrors
%                              - FilteredVariablesKStepAhead
%                              - FilteredVariablesKStepAheadVariances
%    o M_          [structure]  Definition of the model
%    o options_    [structure]  Options; returns options_.first_obs
%    o bayestopt_  [structure]  describing the priors; returns fields like bayestopt_.smoother_var_list from the smoother
%    o Smoothed_variables_declaration_order_deviation_form
%                           Smoothed variables from the Kalman smoother in
%                           order of declaration of variables (M_.endo_names)
%                           in deviations from their respective mean, i.e.
%                           without trend and constant part (used for shock_decomposition)
%    o initial_date [dseries]   initial period, used for shock_decomposition
%
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function use persistent variables for the dataset and the description of the missing observations. Consequently, if this function
%     is called more than once (by changing the value of parameters) the sample *must not* change.

% Copyright © 2010-2020 Dynare Team
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

% store qz_criterium
qz_criterium_old=options_.qz_criterium;

if ischar(parameters) && strcmp(parameters,'calibration')
    options_.smoother=1;
end

[dataset_,dataset_info, ~, ~, M_, options_, oo_, estim_params_,bayestopt_] = dynare_estimation_init(var_list, M_.fname, [], M_, options_, oo_, estim_params_, bayestopt_);
initial_date=dataset_.firstdate; 

% set the qz_criterium
options_=select_qz_criterium_value(options_);

if nargin==0
    parameters = 'posterior_mode';
end

if ischar(parameters)
    switch parameters
      case 'posterior_mode'
        parameters = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
      case 'posterior_mean'
        parameters = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
      case 'posterior_median'
        parameters = get_posterior_parameters('median',M_,estim_params_,oo_,options_);
      case 'mle_mode'
        parameters = get_posterior_parameters('mode',M_,estim_params_,oo_,options_,'mle_');
      case 'prior_mode'
        parameters = bayestopt_.p5(:);
      case 'prior_mean'
        parameters = bayestopt_.p1;
      case 'calibration'
        if isempty(oo_.dr)
            error('You must run ''stoch_simul'' first.');
        end
        parameters = [];
      otherwise
        disp('evaluate_smoother:: If the input argument is a string, then it has to be equal to:')
        disp('                     ''posterior_mode'', ')
        disp('                     ''posterior_mean'', ')
        disp('                     ''posterior_median'', ')
        disp('                     ''prior_mode'' or')
        disp('                     ''prior_mean''.')
        disp('                     ''calibration''.')
        error
    end
end

if options_.occbin.smoother.status
    if options_.occbin.smoother.inversion_filter
        [~, info, ~, ~, ~, ~, ~, ~, ~, ~, oo_, atT, innov] = occbin.IVF_posterior(parameters,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,prior_bounds(bayestopt_,options_.prior_trunc),oo_);
        if ismember(info(1),[303,304,306])
            oo_.occbin.smoother.error_flag=1;
        else
            oo_.occbin.smoother.error_flag=0;
            updated_variables = atT*nan;
            measurement_error=[];
            ys = oo_.dr.ys;
            trend_coeff = zeros(length(options_.varobs_id),1);
            bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
        end
    else
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_] = ...
            occbin.DSGE_smoother(parameters,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info);
    end
else
    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_] = ...
        DsgeSmoother(parameters,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_);
end
if ~(options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter)
    if ~options_.occbin.smoother.status || (options_.occbin.smoother.status && oo_.occbin.smoother.error_flag==0)
        [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
    end
else
    if ~oo_.occbin.smoother.error_flag
        options_nk=options_.nk;
        options_.nk=[]; %unset options_.nk and reset it later
        [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff);
        options_.nk=options_nk;
    else
        fprintf('\nIVF: smoother did not succeed. No results will be written to oo_.\n')
    end
end
if nargout>4
    Smoothed_variables_declaration_order_deviation_form=atT(oo_.dr.inv_order_var(bayestopt_.smoother_var_list),:);
end

%reset qz_criterium
options_.qz_criterium=qz_criterium_old;
oo_.gui.ran_calib_smoother = true;
