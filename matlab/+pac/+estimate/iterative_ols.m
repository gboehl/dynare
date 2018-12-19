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

% Copyright (C) 2018 Dynare Team
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

global M_ oo_

debug = false; % If true, prints the value of SSR for the initial guess.

[pacmodl, ~, rhs, ~, ~, ~, ~, ~, ~, ~, ipnames_, params, data] = ...
    pac.estimate.init(M_, oo_, eqname, params, data, range);

% Set initial condition.
params0 = cell2mat(struct2cell(params));

% Set flag for models with non optimizing agents.
is_non_optimizing_agents = isfield(M_.pac.(pacmodl), 'non_optimizing_behaviour');

if is_non_optimizing_agents
    non_optimizing_behaviour = M_.pac.(pacmodl).non_optimizing_behaviour;
    non_optimizing_behaviour_params = NaN(length(non_optimizing_behaviour.params), 1);
    noparams = isnan(non_optimizing_behaviour.params);
    if ~all(noparams)
        % Find estimated non optimizing behaviour parameters (if any.
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
    estimate_non_optimizing_agents_share = ismember(M_.param_names(M_.pac.(pacmodl).share_of_optimizing_agents_index), fieldnames(params));
    if ~estimate_non_optimizing_agents_share
        share_of_optimizing_agents = M_.params(M_.pac.(pacmodl).share_of_optimizing_agents_index);
        if share_of_optimizing_agents>1 || share_of_optimizing_agents<0
            error('The share of optimizing agents shoud be in (0,1).')
        end
    end
else
    share_of_optimizing_agents = 1.0;
    estimate_non_optimizing_agents_share = false;
end

% Build PAC expectation matrix expression.
dataForPACExpectation0 = dseries();
listofvariables0 = {};
if ~isempty(M_.pac.(pacmodl).h0_param_indices)
    for i=1:length(M_.pac.(pacmodl).h0_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|\\w*\\(-1\\)))\\*%s', M_.param_names{M_.pac.(pacmodl).h0_param_indices(i)}), 'names');
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*\\(-1\\))|(\\w*)))', M_.param_names{M_.pac.(pacmodl).h0_param_indices(i)}), 'names');
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
if ~isempty(M_.pac.(pacmodl).h1_param_indices)
    for i=1:length(M_.pac.(pacmodl).h1_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|(\\w*\\(-1\\))))\\*%s', M_.param_names{M_.pac.(pacmodl).h1_param_indices(i)}), 'names');
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*\\(-1\\))|(\\w*)))', M_.param_names{M_.pac.(pacmodl).h1_param_indices(i)}), 'names');
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
        variable = M_.endo_names{non_optimizing_behaviour.vars(i)};
        if non_optimizing_behaviour.lags(i)
            dataForNonOptimizingBehaviour = [dataForNonOptimizingBehaviour, data{variable}.lag(non_optimizing_behaviour.lags(i))];
        else
            dataForNonOptimizingBehaviour = [dataForNonOptimizingBehaviour, data{variable}];
        end
    end
else
    dataForNonOptimizingBehaviour = dseries();
end

% Reorder ec.vars locally if necessary. Second variable must be the
% endogenous variable, while the first must be the associated trend.
if M_.pac.(pacmodl).ec.isendo(2)
    ecvars = M_.pac.(pacmodl).ec.vars;
else
    ecvars = flip(M_.pac.(pacmodl).ec.vars);
end

% Build matrix for EC and AR terms.
DataForOLS = dseries();

% Error correction term is trend minus the level of the endogenous variable.
DataForOLS{'ec-term'} = data{M_.endo_names{ecvars(1)}}.lag(1)-data{M_.endo_names{ecvars(2)}}.lag(1);
listofvariables3 = {'ec-term'};
xparm = { M_.param_names(M_.pac.(pacmodl).ec.params(1))};
for i = 1:length(M_.pac.(pacmodl).ar.params)
    if islagof(M_.pac.(pacmodl).ar.vars(i), M_.pac.(pacmodl).lhs_var)
        DataForOLS = [DataForOLS, data{M_.endo_names{M_.pac.(pacmodl).ar.vars(i)}}];
        listofvariables3{i+1} = M_.endo_names{M_.pac.(pacmodl).ar.vars(i)};
        xparm{i+1} = M_.param_names(M_.pac.(pacmodl).ar.params(i));
    end
end

XDATA = DataForOLS{listofvariables3{:}}(range).data;

if estimate_non_optimizing_agents_share
    [~, ecm_params_id] = setdiff(ipnames_, M_.pac.(pacmodl).share_of_optimizing_agents_index);
    [~, share_param_id] = setdiff(1:length(ipnames_), ecm_params_id);
    params0_ = params0(ecm_params_id);
    share_of_optimizing_agents = params0(share_param_id);
    if share_of_optimizing_agents>1 || share_of_optimizing_agents<0
        error('Initial value for the share of optimizing agents shoud be in (0,1).')
    end
else
    params0_ = params0;
end

% Update the vector of parameters.
M_.params(ipnames_) = params0;

% Update the reduced form PAC expectation parameters and compute the expectations.
[PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, M_, oo_);

if debug
    YDATA = data{M_.endo_names{M_.pac.(pacmodl).lhs_var}}(range).data;
    if is_non_optimizing_agents
        YDATA = YDATA-share_of_optimizing_agents*PacExpectations;
        YDATA = YDATA-(1-share_of_optimizing_agents)*(dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params);
    else
        YDATA = YDATA-PacExpectations;
    end
    r = YDATA-XDATA*(params0_*share_of_optimizing_agents);
    ssr = r'*r;
    fprintf('\nInitial value of the objective (SSR) is %s.\n\n', num2str(ssr));
end

noconvergence = true;
counter = 0;

while noconvergence
    counter = counter+1;
    % Set vector for left handside variable
    YDATA = data{M_.endo_names{M_.pac.(pacmodl).lhs_var}}(range).data;
    if is_non_optimizing_agents
        YDATA = YDATA-share_of_optimizing_agents*PacExpectations;
        YDATA = YDATA-(1-share_of_optimizing_agents)*(dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params);
    else
        YDATA = YDATA-PacExpectations;
    end
    % Run OLS to estimate PAC parameters (autoregressive parameters and error correction parameter).
    params1_ = (XDATA\YDATA)/share_of_optimizing_agents;
    % Compute residuals and sum of squareed residuals.
    r = YDATA-XDATA*(params1_*share_of_optimizing_agents);
    ssr = r'*r;
    % Update convergence dummy variable and display iteration info.
    noconvergence = max(abs(params0_-params1_))>1e-6;
    fprintf('Iter. %u,\t crit: %10.5f\t ssr: %10.8f\n', counter, max(abs(params0_-params1_)), ssr)
    % Updatevector of estimated parameters.
    params0_ = params1_;
    % Update the value of the share of non optimizing agents (if estimated)
    if estimate_non_optimizing_agents_share
        % First update the parameters and compute the PAC expectation reduced form parameters.
        M_.params(ipnames_(ecm_params_id)) = params0_;
        [PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, M_, oo_);
        % Set vector for left handside variable.
        YDATA = data{M_.endo_names{M_.pac.(pacmodl).lhs_var}}(range).data;
        YDATA = YDATA-dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params;
        % Set vector for regressor.
        ZDATA = XDATA*params0_+PacExpectations-dataForNonOptimizingBehaviour(range).data*non_optimizing_behaviour_params;
        % Update the (estimated) share of optimizing agents by running OLS
        share_of_optimizing_agents = (ZDATA\YDATA);
        % Force the share of optimizing agents to be in [0,1].
        share_of_optimizing_agents = max(min(share_of_optimizing_agents, 1.0), 0.0);
        M_.params(ipnames_(share_param_id)) = share_of_optimizing_agents;
    else
        M_.params(ipnames_) = params0_;
        [PacExpectations, M_] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, M_, oo_);
    end
end



function [PacExpectations, Model] = UpdatePacExpectationsData(dataPAC0, dataPAC1, data, range, pacmodl, Model, Output)
    % Update PAC reduced parameters.
    Model = pac.update.parameters(pacmodl, Model, Output);
    % Compute PAC expectation.
    if isempty(dataPAC0)
        PacExpectations = 0;
    else
        PacExpectations = dataPAC0*Model.params(Model.pac.pacman.h0_param_indices);
    end
    if ~isempty(dataPAC1)
        PacExpectations = PacExpectations+dataPAC1*Model.params(Model.pac.pacman.h1_param_indices);
    end
    % Add correction for growth neutrality if required.
    correction = 0;
    if isfield(Model.pac.(pacmodl), 'growth_type')
        switch Model.pac.(pacmodl).growth_type
          case 'parameter'
            correction = Model.params(Model.pac.(pacmodl).growth_index)*Model.params(Model.pac.(pacmodl).growth_neutrality_param_index);
          case 'exogenous'
            GrowthVariable = data{Model.exo_names{Model.pac.(pacmodl).growth_index}};
            GrowthVariable = GrowthVariable(range).data;
            correction = GrowthVariable*Model.params(Model.pac.(pacmodl).growth_neutrality_param_index);
          case 'endogenous'
            GrowthVariable = data{Model.endo_names{Model.pac.(pacmodl).growth_index}};
            GrowthVariable = GrowthVariable(range).data;
            correction = GrowthVariable*Model.params(Model.pac.(pacmodl).growth_neutrality_param_index);
          otherwise
            error('Not yet implemented.')
        end
    end
    PacExpectations = PacExpectations+correction;