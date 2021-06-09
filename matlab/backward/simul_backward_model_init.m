function [initialconditions, samplesize, innovations, DynareOptions, DynareModel, DynareOutput, endonames, exonames, nx, ny1, iy1, jdx, model_dynamic, y] = ...
    simul_backward_model_init(initialconditions, samplesize, DynareOptions, DynareModel, DynareOutput, innovations)

% Initialization of the routines simulating backward models.

% Copyright (C) 2017-2019 Dynare Team
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

% Test if the model is backward.
if DynareModel.maximum_lead
    error('simul_backward_nonlinear_model:: The specified model is not backward looking!')
end

% Test if the first argument is a dseries object.
if ~(isdseries(initialconditions) || isempty(initialconditions))
    error('First input argument must be a dseries object or an empty array!')
end

% If initialconditions is empty instantiates a dseries object with the informations available in DynareModel.endo_histval.
if isempty(initialconditions)
    yinitdata = zeros(DynareModel.orig_endo_nbr, DynareModel.orig_maximum_lag);
    yinitdata(:,1) = DynareModel.endo_histval(1:DynareModel.orig_endo_nbr);
    xinitdata = zeros(DynareModel.exo_nbr, DynareModel.orig_maximum_lag);
    if DynareModel.orig_maximum_endo_lag>1
        for i = 1:length(DynareModel.aux_vars)
            if DynareModel.aux_vars(i).type==1
                yinitdata(DynareModel.aux_vars(i).orig_index, abs(DynareModel.aux_vars(i).orig_lead_lag)+1) = ...
                    DynareModel.endo_histval(DynareModel.orig_endo_nbr+i);
            end
        end
        yinitdata = flip(yinitdata, 2);
    end
    if DynareModel.orig_maximum_exo_lag>0
        for i = 1:length(DynareModel.aux_vars)
            if DynareModel.aux_vars(i).type==3
                xinitdata(DynareModel.aux_vars(i).orig_index, abs(DynareModel.aux_vars(i).orig_lead_lag)+1) = ...
                    DynareModel.endo_histval(DynareModel.orig_endo_nbr+i);
            end
        end
        xinitdata = flip(xinitdata, 2);
    end
    initialconditions = dseries([transpose(yinitdata) transpose(xinitdata)], '1Y', ...
                                vertcat(DynareModel.endo_names(1:DynareModel.orig_endo_nbr), DynareModel.exo_names));
end

[initialconditions, info] = checkdatabase(initialconditions, DynareModel, false, true);

% Test if the first argument contains all the lagged endogenous variables
endonames = DynareModel.endo_names;
missingendogenousvariables = setdiff(endonames, initialconditions.name);
endolags = get_lags_on_endogenous_variables(DynareModel);
endolags_ = endolags(find(endolags));
endowithlagnames = endonames(find(endolags));
if ~isempty(missingendogenousvariables)
    missingendogenousvariables = setdiff(endowithlagnames, initialconditions.name);
    missingendogenouslaggedvariables = intersect(endowithlagnames, missingendogenousvariables);
    if ~isempty(missingendogenouslaggedvariables)
        disp('You have to initialize the following endogenous variables:')
        msg = sprintf('%s\n', missingendogenouslaggedvariables{1:end-1});
        msg = sprintf('%s%s', msg, missingendogenouslaggedvariables{end});
        disp(msg)
        skipline()
        error('Please fix the dseries object used for setting the initial conditions!')
    end
end

% Test if we have enough periods in the database.
maxlag = abs(min(endolags));
if maxlag>initialconditions.nobs
    error('The dseries object provided as first input argument should at least have %s periods!', num2str(maxlag))
end
missinginitialcondition = false;
for i = 1:length(endowithlagnames)
    lags = abs(endolags_(i));
    variable = initialconditions{endowithlagnames{i}};
    nanvalues = isnan(variable.data);
    if any(nanvalues(end-(lags-1):end))
        missinginitialcondition = true;
        for j=variable.nobs:-1:variable.nobs-(lags-1)
            if isnan(variable.data(j))
                dprintf('Variable %s should not have a NaN value in period %s.', endowithlagnames{i}, date2string(variable.dates(j)))
            end
        end
    end
end
if missinginitialcondition
    skipline()
    error('Please fix the dseries object used for setting the initial conditions!')
end

% If the model has lags on the exogenous variables, test if we have corresponding initial conditions.
exonames = DynareModel.exo_names;
missingexogenousvariables = setdiff(exonames, initialconditions.name);
exolags = get_lags_on_exogenous_variables(DynareModel);
exolags_ = exolags(find(exolags));
exowithlagnames = exonames(find(exolags));
if ~isempty(missingexogenousvariables)
    missingexogenousvariables = setdiff(exowithlagnames, initialconditions.name);
    missingexogenouslaggedvariables = intersect(exowithlagnames, missingexogenousvariables);
    if ~isempty(missingexogenouslaggedvariables)
        disp('You have to initialize the following exogenous variables:')
        msg = sprintf('%s\n', missingexogenouslaggedvariables{1:end-1});
        msg = sprintf('%s%s', msg, missingexogenouslaggedvariables{end});
        disp(msg)
        skipline()
        error('Please fix the dseries object used for setting the initial conditions!')
    end
end

% Test if we have enough periods in the database.
maxlag = abs(min(exolags));
if maxlag>initialconditions.nobs
    error('The dseries object provided as first input argument should at least have %s periods!', num2str(maxlag))
end
missinginitialcondition = false;
for i = 1:length(exowithlagnames)
    lags = abs(exolags_(i));
    variable = initialconditions{exowithlagnames{i}};
    nanvalues = isnan(variable.data);
    if any(nanvalues(end-(lags-1):end))
        missinginitialcondition = true;
        for j=variable.nobs:-1:variable.nobs-(lags-1)
            if isnan(variable.data(j))
                dprintf('Variable %s should not have a NaN value in period %s.', exowithlagnames{i}, date2string(variable.dates(j)))
            end
        end
    end
end
if missinginitialcondition
    skipline()
    error('Please fix the dseries object used for setting the initial conditions!')
end

if nargin<6 || isempty(innovations)
    % Set the covariance matrix of the structural innovations.
    variances = diag(DynareModel.Sigma_e);
    number_of_shocks = length(DynareModel.Sigma_e);
    positive_var_indx = find(variances>0);
    effective_number_of_shocks = length(positive_var_indx);
    covariance_matrix = DynareModel.Sigma_e(positive_var_indx,positive_var_indx);
    covariance_matrix_upper_cholesky = chol(covariance_matrix);
    % Set seed to its default state.
    if DynareOptions.bnlms.set_dynare_seed_to_default
        set_dynare_seed('default');
    end
    % Simulate structural innovations.
    switch DynareOptions.bnlms.innovation_distribution
      case 'gaussian'
        DynareOutput.bnlms.shocks = randn(samplesize,effective_number_of_shocks)*covariance_matrix_upper_cholesky;
      otherwise
        error(['simul_backward_nonlinear_model:: ' DynareOption.bnlms.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
    % Put the simulated innovations in DynareOutput.exo_simul.
    DynareOutput.exo_simul = zeros(samplesize,number_of_shocks);
    DynareOutput.exo_simul(:,positive_var_indx) = DynareOutput.bnlms.shocks;
    innovations = DynareOutput.exo_simul;
else
    DynareOutput.exo_simul = innovations; % innovations
end

% Initialization of the returned simulations.
DynareOutput.endo_simul = NaN(DynareModel.endo_nbr, samplesize+initialconditions.nobs);
for i=1:length(endonames)
    if ismember(endonames{i}, initialconditions.name)
        DynareOutput.endo_simul(i,1:initialconditions.nobs) = transpose(initialconditions{endonames{i}}.data);
    end
end

% Initialization of the array for the exogenous variables.
DynareOutput.exo_simul = [NaN(initialconditions.nobs, DynareModel.exo_nbr); DynareOutput.exo_simul ];
for i=1:length(exonames)
    if ismember(exonames{i}, initialconditions.name)
        DynareOutput.exo_simul(1:initialconditions.nobs, i) = initialconditions{exonames{i}}.data;
    end
end


if nargout>8
    nx = size(DynareOutput.exo_simul,2);
    ny0 = nnz(DynareModel.lead_lag_incidence(2,:));
    ny1 = nnz(DynareModel.lead_lag_incidence(1,:));
    iy1 = find(DynareModel.lead_lag_incidence(1,:)>0);
    idx = 1:DynareModel.endo_nbr;
    jdx = idx+ny1;
    % Get the name of the dynamic model routine.
    model_dynamic = str2func([DynareModel.fname,'.dynamic']);
    % initialization of vector y.
    y = NaN(length(idx)+ny1,1);
end