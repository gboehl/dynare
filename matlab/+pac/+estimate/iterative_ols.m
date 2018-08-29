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

[pacmodl, lhs, rhs, pnames, enames, xnames, pid, eid, xid, pnames_, ipnames_, params, data, islaggedvariables] = ...
    pac.estimate.init(M_, oo_, eqname, params, data, range);

% Set initial condition.
params0 = cell2mat(struct2cell(params));

% Build PAC expectation matrix expression.
dataForPACExpectation0 = dseries();
listofvariables0 = {};
if ~isempty(M_.pac.(pacmodl).h0_param_indices)
    for i=1:length(M_.pac.(pacmodl).h0_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|\\w*\\(-1\\)))\\*%s', M_.param_names{M_.pac.(pacmodl).h0_param_indices(i)}), 'names'); 
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*)|\\w*\\(-1\\)))', M_.param_names{M_.pac.(pacmodl).h0_param_indices(i)}), 'names'); 
        end
        if isempty(strfind(match.var, '(-1)'))
            listofvariables0{i} = match.var;
            dataForPACExpectation0 = [dataForPACExpectation0, data{listofvariables0{i}}];
        else
            listofvariables0{i} = match.var(1:end-4);
            dataForPACExpectation0 = [dataForPACExpectation0, data{match.var(1:end-4)}.lag(1)];
        end
    end
end
dataForPACExpectation1 = dseries();
listofvariables1 = {};
if ~isempty(M_.pac.(pacmodl).h1_param_indices)
    for i=1:length(M_.pac.(pacmodl).h1_param_indices)
        match = regexp(rhs, sprintf('(?<var>((\\w*)|\\w*\\(-1\\)))\\*%s', M_.param_names{M_.pac.(pacmodl).h1_param_indices(i)}), 'names'); 
        if isempty(match)
            match = regexp(rhs, sprintf('%s\\*(?<var>((\\w*)|\\w*\\(-1\\)))', M_.param_names{M_.pac.(pacmodl).h1_param_indices(i)}), 'names'); 
        end
        if isempty(strfind(match.var, '(-1)'))
            listofvariables1{i} = match.var;
            dataForPACExpectation1 = [dataForPACExpectation1, data{listofvariables1{i}}];
        else
            listofvariables1{i} = match.var(1:end-4);
            dataForPACExpectation1 = [dataForPACExpectation1, data{match.var(1:end-4)}.lag(1)];
        end
    end
end

% Set correction for growth neutrality
correction = 0;
if isfield(M_.pac.(pacmodl), 'growth_type')
    switch M_.pac.(pacmodl).growth_type
      case 'parameter'
        correction = M_.params(M_.pac.(pacmodl).growth_index)*M_.params(M_.pac.(pacmodl).growth_neutrality_param_index);
      case 'exogenous'
        ExogenousGrowthVariable = ...
            data{M_.exo_names{M_.pac.(pacmodl).growth_index}};
        ExogenousGrowthVariable = ExogenousGrowthVariable(range).data;
        correction = ExogenousGrowthVariable*M_.params(M_.pac.(pacmodl).growth_neutrality_param_index);
      otherwise
        error('Not yet implemented.')
    end
end

% Build matrix for EC and AR terms.
DataForOLS = dseries();
DataForOLS{'ec-term'} = data{M_.endo_names{M_.pac.(pacmodl).ec.vars(1)}}.lag(1)-data{M_.endo_names{M_.pac.(pacmodl).ec.vars(2)}}.lag(1);
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

params0

% Iterative OLS
noconvergence = true;
counter = 0;
while noconvergence
    counter = counter+1;
    % Update the PAC parameter,s.
    for i=1:length(ipnames_)
        M_.params(ipnames_(i)) = params0(i);
    end
    M_ = pac.update.parameters(pacmodl, M_, oo_);
    % Compute lhs conditional on params0
    PACExpectations = 0;
    if ~isempty(listofvariables0)
        PACExpectations = dataForPACExpectation0{listofvariables0{:}}(range).data*M_.params(M_.pac.pacman.h0_param_indices);
    end
    if ~isempty(listofvariables1)
        PACExpectations = dataForPACExpectation1{listofvariables1{:}}(range).data*M_.params(M_.pac.pacman.h1_param_indices);
    end
    YDATA = data{M_.endo_names{M_.pac.(pacmodl).lhs_var}}(range).data-correction-PACExpectations;
    % Do OLS
    params1 = XDATA\YDATA;
    noconvergence = max(abs(params0-params1))>1e-6;
    disp(sprintf('Iter. %u,\t crit: %10.5f', counter, max(abs(params0-params1))))
    params0 = params1;
end


% Update M_.params
for i=1:length(params0)
    M_.params(ipnames_(i)) = params0(i);
end

M_ = pac.update.parameters(pacmodl, M_, oo_);