function iterative_ols(eqname, params, data, range)

% Estimates the parameters of a PAC equation by Iterative Ordinary Least Squares.
%
% INPUTS
% - eqname       [string]    Name of the pac equation.
% - params       [struct]    Describes the parameters to be estimated.
% - data         [dseries]   Database for the estimation
% - range        [dates]     Range of dates for the estimation.
%
% OUTPUTS
% - none
%
% REMARKS
% [1] The estimation results are printed in the command line window, and the
%     parameters are updated accordingly in M_.params.
% [2] The second input is a structure. Each fieldname corresponds to the
%     name of an estimated parameter, the value of the field is the initial
%     guess used for the estimation (by NLS).
% [3] The third input is a dseries object which must at least contain all
%     the variables appearing in the estimated equation. The residual of the
%     equation must have NaN values in the object.
% [4] It is assumed that the residual is additive.

% Copyright (C) 2018-2021 Dynare Team
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

global M_ oo_ options_

[pacmodl, ~, rhs, ~, ~, ~, rname, ~, ~, ~, ~, ipnames_, params, data, ~, eqtag] = ...
    pac.estimate.init(M_, oo_, eqname, params, data, range);

% Set initial condition.
params0 = cell2mat(struct2cell(params));

% Set flag for models with non optimizing agents.
is_non_optimizing_agents = isfield(M_.pac.(pacmodl).equations.(eqtag), 'non_optimizing_behaviour');

% Set flag for models with exogenous variables (outside of non optimizing agents part)
if isfield(M_.pac.(pacmodl).equations.(eqtag), 'additive')
    is_exogenous_variables = length(M_.pac.(pacmodl).equations.(eqtag).additive.vars)>1;
else
    is_exogenous_variables = false;
end

% Set flag for models with exogenous variables (in the optimizing agents part)
if isfield(M_.pac.(pacmodl).equations.(eqtag), 'optim_additive')
    is_optim_exogenous_variables = length(M_.pac.(pacmodl).equations.(eqtag).optim_additive.vars)>0;
else
    is_optim_exogenous_variables = false;
end

if is_non_optimizing_agents
    non_optimizing_behaviour = M_.pac.(pacmodl).equations.(eqtag).non_optimizing_behaviour;
    non_optimizing_behaviour_params = NaN(length(non_optimizing_behaviour.params), 1);
    noparams = isnan(non_optimizing_behaviour.params);
    if ~all(noparams)
        % Find estimated non optimizing behaviour parameters (if any).
        non_optimizing_behaviour_estimated_params = ismember(M_.param_names(non_optimizing_behaviour.params), fieldnames(params));
        if any(non_optimizing_behaviour_estimated_params)
            error('The estimation of non optimizing behaviour parameters is not yet allowed.')
        else
            non_optimizing_behaviour_params(noparams) = 1.0;
            non_optimizing_behaviour_params(~noparams) = M_.params(non_optimizing_behaviour.params(~noparams));
        end
    else
        non_optimizing_behaviour_params(noparams) = 1.0;
    end
    non_optimizing_behaviour_params = non_optimizing_behaviour_params.*transpose(non_optimizing_behaviour.scaling_factor);
    % Set flag for the estimation of the share of non optimizing agents.
    estimate_non_optimizing_agents_share = ismember(M_.param_names(M_.pac.(pacmodl).equations.(eqtag).share_of_optimizing_agents_index), fieldnames(params));
    if ~estimate_non_optimizing_agents_share
        share_of_optimizing_agents = M_.params(M_.pac.(pacmodl).equations.(eqtag).share_of_optimizing_agents_index);
        if share_of_optimizing_agents>1 || share_of_optimizing_agents<0
            error('The share of optimizing agents shoud be in (0,1).')
        end
    end
    share_of_optimizing_agents_index = M_.pac.(pacmodl).equations.(eqtag).share_of_optimizing_agents_index;
else
    share_of_optimizing_agents = 1.0;
    share_of_optimizing_agents_index = [];
    estimate_non_optimizing_agents_share = false;
end

if is_exogenous_variables
    additive = M_.pac.(pacmodl).equations.(eqtag).additive;
    residual_id = find(strcmp(rname, M_.exo_names));
    residual_jd = find(additive.vars==residual_id & ~additive.isendo);
    additive.params(residual_jd) = [];
    additive.vars(residual_jd) = [];
    additive.isendo(residual_jd) = [];
    additive.lags(residual_jd) = [];
    additive.scaling_factor(residual_jd) = [];
    additive.estimation  = ismember(additive.params, ipnames_);
else
    additive = struct('params', [], 'vars', [], 'isendo', [], 'lags', [], 'scaling_factor', [], 'estimation', []);
end

if is_optim_exogenous_variables
    optim_additive = M_.pac.(pacmodl).equations.(eqtag).optim_additive;
    optim_additive.estimation  = ismember(optim_additive.params, ipnames_);
else
    optim_additive = struct('params', [], 'vars', [], 'isendo', [], 'lags', [], 'scaling_factor', [], 'estimation', []);
end


% Build PAC expectation matrix expression.
dataForPACExpectation0 = dseries();
listofvariables0 = {};
if ~isempty(M_.pac.(pacmodl).equations.(eqtag).h0_param_indices)
    for i=1:length(M_.pac.(pacmodl).equations.(eqtag).h0_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|\\w*\\(-1\\)))\\*%s', M_.param_names{M_.pac.(pacmodl).equations.(eqtag).h0_param_indices(i)}), 'names');
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*\\(-1\\))|(\\w*)))', M_.param_names{M_.pac.(pacmodl).equations.(eqtag).h0_param_indices(i)}), 'names');
        end
        if isempty(strfind(match.var, '(-1)'))
            listofvariables0{i} = match.var;
            dataForPACExpectation0 = [dataForPACExpectation0, data{listofvariables0{i}}];
        else
            listofvariables0{i} = match.var(1:end-4);
            dataForPACExpectation0 = [dataForPACExpectation0, data{match.var(1:end-4)}.lag(1)];
        end
    end
    dataPAC0 = dataForPACExpectation0{listofvariables0{:}}(range).data;
else
    dataPAC0 = [];
end

dataForPACExpectation1 = dseries();
listofvariables1 = {};
if ~isempty(M_.pac.(pacmodl).equations.(eqtag).h1_param_indices)
    for i=1:length(M_.pac.(pacmodl).equations.(eqtag).h1_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|(\\w*\\(-1\\))))\\*%s', M_.param_names{M_.pac.(pacmodl).equations.(eqtag).h1_param_indices(i)}), 'names');
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*\\(-1\\))|(\\w*)))', M_.param_names{M_.pac.(pacmodl).equations.(eqtag).h1_param_indices(i)}), 'names');
        end
        if isempty(strfind(match.var, '(-1)'))
            listofvariables1{i} = match.var;
            dataForPACExpectation1 = [dataForPACExpectation1, data{listofvariables1{i}}];
        else
            listofvariables1{i} = match.var(1:end-4);
            dataForPACExpectation1 = [dataForPACExpectation1, data{match.var(1:end-4)}.lag(1)];
        end
    end
    dataPAC1 = dataForPACExpectation1{listofvariables1{:}}(range).data;
else
    dataPAC1 = [];
end

% Build data for non optimizing behaviour
if is_non_optimizing_agents
    dataForNonOptimizingBehaviour = dseries();
    for i=1:length(non_optimizing_behaviour.vars)
        if non_optimizing_behaviour.isendo(i)
            variable = M_.endo_names{non_optimizing_behaviour.vars(i)};
        else
            variable = M_.exo_names{non_optimizing_behaviour.vars(i)};
        end
        if non_optimizing_behaviour.lags(i)
            dataForNonOptimizingBehaviour = [dataForNonOptimizingBehaviour, data{variable}.lag(non_optimizing_behaviour.lags(i))];
        else
            dataForNonOptimizingBehaviour = [dataForNonOptimizingBehaviour, data{variable}];
        end
    end
else
    dataForNonOptimizingBehaviour = dseries();
end

% Build data for exogenous variables (out of non optimizing behaviour term).
if is_exogenous_variables
    listofvariables2 = {}; j = 0;
    dataForExogenousVariables = dseries();  % Estimated parameters
    dataForExogenousVariables_ = 0;         % Calibrated parameters
    is_any_calibrated_parameter_x = false;
    is_any_estimated_parameter_x = false;
    for i=1:length(additive.vars)
        if additive.isendo(i)
            variable = M_.endo_names{additive.vars(i)};
        else
            variable = M_.exo_names{additive.vars(i)};
        end
        if additive.estimation(i)
            j = j+1;
            is_any_estimated_parameter_x = true;
            listofvariables2{j} = variable;
            dataForExogenousVariables = [dataForExogenousVariables, additive.scaling_factor(i)*data{variable}.lag(-additive.lags(i))];
        else
            is_any_calibrated_parameter_x = true;
            tmp = data{variable}.lag(-additive.lags(i)).data;
            if ~isnan(additive.params(i))
                tmp = M_.params(additive.params(i))*tmp;
            end
            tmp = additive.scaling_factor(i)*tmp;
            dataForExogenousVariables_ = dataForExogenousVariables_+tmp;
        end
    end
    if is_any_calibrated_parameter_x
        dataForExogenousVariables_ = dseries(dataForExogenousVariables_, data.dates(1), 'exogenous_variables_associated_with_calibrated_parameters');
    end
else
    listofvariables2 = {};
    dataForExogenousVariables = dseries();
    dataForExogenousVariables_ = dseries();
end

% Build data for exogenous variables (in the optimizing behaviour term).
if is_optim_exogenous_variables
    listofvariables4 = {}; j = 0;
    dataForOptimExogenousVariables = dseries();  % Estimated parameters
    dataForOptimExogenousVariables_ = 0;         % Calibrated parameters
    is_any_calibrated_parameter_optim_x = false;
    is_any_estimated_parameter_optim_x = false;
    for i=1:length(optim_additive.vars)
        if optim_additive.isendo(i)
            variable = M_.endo_names{optim_additive.vars(i)};
        else
            variable = M_.exo_names{optim_additive.vars(i)};
        end
        if optim_additive.estimation(i)
            j = j+1;
            is_any_estimated_parameter_optim_x = true;
            listofvariables4{j} = variable;
            dataForOptimExogenousVariables = [dataForOptimExogenousVariables, optim_additive.scaling_factor(i)*data{variable}.lag(-optim_additive.lags(i))];
        else
            is_any_calibrated_parameter_optim_x = true;
            tmp = data{variable}.lag(-optim_additive.lags(i)).data;
            if ~isnan(optim_additive.params(i))
                tmp = M_.params(optim_additive.params(i))*tmp;
            end
            tmp = optim_additive.scaling_factor(i)*tmp;
            dataForOptimExogenousVariables_ = dataForOptimExogenousVariables_+tmp;
        end
    end
    if is_any_calibrated_parameter_optim_x
        dataForOptimExogenousVariables_ = dseries(dataForOptimExogenousVariables_, data.dates(1), 'exogenous_variables_associated_with_calibrated_parameters');
    end
else
    listofvariables4 = {};
    dataForOptimExogenousVariables = dseries();
    dataForOptimExogenousVariables_ = dseries();
end

% Reorder ec.vars locally if necessary. Second variable must be the
% endogenous variable, while the first must be the associated trend.
if M_.pac.(pacmodl).equations.(eqtag).ec.isendo(2)
    ecvars = M_.pac.(pacmodl).equations.(eqtag).ec.vars;
else
    ecvars = flip(M_.pac.(pacmodl).equations.(eqtag).ec.vars);
end

%% Build matrix for EC and AR terms.
DataForOLS = dseries();

% Error correction term is trend minus the level of the endogenous variable.
DataForOLS{'ec-term'} = data{M_.endo_names{ecvars(1)}}.lag(1)-data{M_.endo_names{ecvars(2)}}.lag(1);
listofvariables3 = {'ec-term'};
xparm = { M_.param_names(M_.pac.(pacmodl).equations.(eqtag).ec.params(1))};
for i = 1:length(M_.pac.(pacmodl).equations.(eqtag).ar.params)
    if islagof(M_.pac.(pacmodl).equations.(eqtag).ar.vars(i), M_.pac.(pacmodl).equations.(eqtag).lhs_var)
        DataForOLS = [DataForOLS, data{M_.endo_names{M_.pac.(pacmodl).equations.(eqtag).ar.vars(i)}}];
        listofvariables3{i+1} = M_.endo_names{M_.pac.(pacmodl).equations.(eqtag).ar.vars(i)};
        xparm{i+1} = M_.param_names(M_.pac.(pacmodl).equations.(eqtag).ar.params(i));
    end
end

XDATA = DataForOLS{listofvariables3{:}}(range).data;

if is_optim_exogenous_variables && is_any_estimated_parameter_optim_x
    XDATA = [XDATA, dataForOptimExogenousVariables{listofvariables4{:}}(range).data];
end

if is_exogenous_variables && is_any_estimated_parameter_x
    XDATA = [XDATA, dataForExogenousVariables{listofvariables2{:}}(range).data];
end

% Get index in params0 for share of optimizing agents parameter (if
% not estimated, params_id_0 is empty).
if is_non_optimizing_agents
    params_id_0 = find(ipnames_==share_of_optimizing_agents_index);
else
    params_id_0 = [];
end

% Get indices in params0 for EC and AR parameters
[~, params_id_1] = setdiff(ipnames_, [share_of_optimizing_agents_index, optim_additive.params, additive.params, ]);

% Get indices in params0 for EC and AR parameters plus parameters related to exogenous variables in the optimal part.
[~, params_id_5] = setdiff(ipnames_, [share_of_optimizing_agents_index, additive.params]);

% Get indices in params0 for other parameters (optimizing agents share plus parameters related to exogenous variables).
[~, params_id_2] = setdiff(1:length(ipnames_), params_id_1);
if isoctave && octave_ver_less_than('6')
    % Under Octave < 6, setdiff() behaves as with the 'legacy' option (i.e. pre-R2012b behaviour)
    params_id_2 = params_id_2';
end

% Get indices in params0 for the parameters associated to the exogenous variables.
params_id_3 = setdiff(params_id_2, params_id_0);
[~, params_id_3_o] = ismember(optim_additive.params, ipnames_);
params_id_3_no = setdiff(params_id_3, params_id_3_o);

% Get indices in params0 for EC/AR parameters and parameters associated to the exogenous variables (if any).
params_id_4 = [params_id_1; params_id_3];

% Get values for EC-AR parameters plus the parameters associated to the exogenous variables (if any).
params0_ = params0([params_id_1; params_id_3]);

% Get value of the share of optimizing agents.
if estimate_non_optimizing_agents_share
    share_of_optimizing_agents = params0(params_id_0);
end

% Check that the share is in (0,1)
if share_of_optimizing_agents>1 || share_of_optimizing_agents<0
    error('Initial value for the share of optimizing agents shoud be in (0,1).')
end

% Update the vector of parameters.
M_.params(ipnames_) = params0;

% Update the reduced form PAC expectation parameters and compute the expectations.
[PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, eqtag, M_, oo_);

noconvergence = true;
counter = 0;

while noconvergence
    counter = counter+1;
    % Set vector for left handside variable
    YDATA = data{M_.endo_names{M_.pac.(pacmodl).equations.(eqtag).lhs_var}}(range).data;
    if is_non_optimizing_agents
        YDATA = YDATA-share_of_optimizing_agents*PacExpectations;
        YDATA = YDATA-(1-share_of_optimizing_agents)*(dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params);
    else
        YDATA = YDATA-PacExpectations;
    end
    if is_exogenous_variables && is_any_calibrated_parameter_x
        YDATA = YDATA-dataForExogenousVariables_(range).data;
    end
    if is_optim_exogenous_variables && is_any_calibrated_parameter_optim_x
        YDATA = YDATA-share_of_optimizing_agents*dataForOptimExogenousVariables_(range).data;
    end
    % Run OLS to estimate PAC parameters (autoregressive parameters and error correction parameter).
    params1_ = (XDATA\YDATA);
    params1_(1:length(params_id_5)) = params1_(1:length(params_id_5))/share_of_optimizing_agents;
    % Compute residuals and sum of squareed residuals.
    r = YDATA-XDATA(:,1:length(params_id_5))*params1_(1:length(params_id_5))*share_of_optimizing_agents;
    if is_optim_exogenous_variables && is_any_estimated_parameter_optim_x
        r = r-XDATA(:,length(params_id_5)+1:end)*params1_(length(params_id_5)+1:end);
    end
    ssr = r'*r;
    % Update convergence dummy variable and display iteration info.
    noconvergence = max(abs(params0_-params1_))>1e-6;
    fprintf('Iter. %u,\t crit: %10.5f\t ssr: %10.8f\n', counter, max(abs(params0_-params1_)), ssr)
    % Updatevector of estimated parameters.
    params0_ = params1_;
    % Update the value of the share of non optimizing agents (if estimated)
    if estimate_non_optimizing_agents_share
        % First update the parameters and compute the PAC expectation reduced form parameters.
        M_.params(ipnames_(params_id_4)) = params0_;
        [PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, eqtag, M_, oo_);
        % Set vector for left handside variable.
        YDATA = data{M_.endo_names{M_.pac.(pacmodl).equations.(eqtag).lhs_var}}(range).data;
        YDATA = YDATA-dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params;
        if is_exogenous_variables
            if is_any_calibrated_parameter_x
                YDATA = YDATA-dataForExogenousVariables_(range).data;
            end
            if is_any_estimated_parameter_x
               YDATA = YDATA-XDATA(:, length(params_id_1):end)*params0_(length(params_id_1):end);
            end
        end
        % Set vector for regressor.
        ZDATA = XDATA(:,1:length(params_id_5))*params0_(1:length(params_id_5))+PacExpectations-dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params;
        if is_optim_exogenous_variables && is_any_calibrated_parameter_optim_x
            ZDATA = ZDATA+dataForOptimExogenousVariables_(range).data;
        end
        % Update the (estimated) share of optimizing agents by running OLS
        share_of_optimizing_agents = (ZDATA\YDATA);
        % Force the share of optimizing agents to be in [0,1].
        share_of_optimizing_agents = max(min(share_of_optimizing_agents, options_.pac.estimation.ols.share_of_optimizing_agents.ub), ...
                                         options_.pac.estimation.ols.share_of_optimizing_agents.lb);
        % Issue an error if the share is nor strictly positive
        if share_of_optimizing_agents<eps()
            error('On iteration %u the share of optimizing agents is found to be zero. Please increase the default lower bound for this parameter.', counter)
        end
        M_.params(ipnames_(params_id_0)) = share_of_optimizing_agents;
    else
        M_.params(ipnames_) = params0_;
        [PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, eqtag, M_, oo_);
    end
end

% Save results
oo_.pac.(pacmodl).equations.(eqtag).ssr = ssr;
oo_.pac.(pacmodl).equations.(eqtag).residual = r;
oo_.pac.(pacmodl).equations.(eqtag).estimator = params0_;
oo_.pac.(pacmodl).equations.(eqtag).parnames = fieldnames(params);
oo_.pac.(pacmodl).equations.(eqtag).covariance = NaN(length(params0_));
oo_.pac.(pacmodl).equations.(eqtag).student = NaN(size(params0_));


function [PacExpectations, Model] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, eqtag, Model, Output)
    % Update PAC reduced parameters.
    Model = pac.update.parameters(pacmodl, Model, Output, false);
    % Compute PAC expectation.
    if isempty(dataPAC0)
        PacExpectations = 0;
    else
        PacExpectations = dataPAC0*Model.params(Model.pac.(pacmodl).equations.(eqtag).h0_param_indices);
    end
    if ~isempty(dataPAC1)
        PacExpectations = PacExpectations+dataPAC1*Model.params(Model.pac.(pacmodl).equations.(eqtag).h1_param_indices);
    end
    % Add correction for growth neutrality if required.
    correction = 0;
    if isfield(Model.pac.(pacmodl), 'growth_linear_comb')
        for iter = 1:numel(Model.pac.(pacmodl).growth_linear_comb)
            GrowthVariable = Model.pac.(pacmodl).growth_linear_comb(iter).constant;
            if Model.pac.(pacmodl).growth_linear_comb(iter).param_id > 0
                GrowthVariable = GrowthVariable*Model.params(Model.pac.(pacmodl).growth_linear_comb(iter).param_id);
            end
            if Model.pac.(pacmodl).growth_linear_comb(iter).exo_id > 0
                GrowthVariable = GrowthVariable*data{Model.exo_names{Model.pac.(pacmodl).growth_linear_comb(iter).exo_id}}.lag(abs(Model.pac.(pacmodl).growth_linear_comb(iter).lag));
                GrowthVariable = GrowthVariable(range).data;
            elseif Model.pac.(pacmodl).growth_linear_comb(iter).endo_id > 0
                GrowthVariable = GrowthVariable*data{Model.endo_names{Model.pac.(pacmodl).growth_linear_comb(iter).endo_id}}.lag(abs(Model.pac.(pacmodl).growth_linear_comb(iter).lag));
                GrowthVariable = GrowthVariable(range).data;
            end
            correction = correction + GrowthVariable;
        end
        correction = correction*Model.params(Model.pac.(pacmodl).growth_neutrality_param_index);
    end
    PacExpectations = PacExpectations+correction;
