function M_ = update_parameters(varexpectationmodelname, M_, oo_)
% function M_ = update_parameters(varexpectationmodelname, M_, oo_)
% Updates the VAR expectation reduced form parameters.
%
% INPUTS
% - varexpectationmodelname       [string]    Name of the pac equation.
% - M_                            [struct]    global structure (model properties)
% - oo_                           [struct]    oo_ global structure (model results)
%
% OUTPUTS
% - M_                            [struct]    M_ global structure (with updated params field)
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2018-2023 Dynare Team
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

% Check that the first input is a row character array.
if ~isrow(varexpectationmodelname)==1 || ~ischar(varexpectationmodelname)
    error('First input argument must be a row character array!')
end

% Check that the model exists.
if ~isfield(M_.var_expectation, varexpectationmodelname)
    error('VAR_EXPECTATION_MODEL %s is not defined!', varexpectationmodelname)
end

% Get the VAR model description
varexpectationmodel = M_.var_expectation.(varexpectationmodelname);

% Get the name of the associated VAR model and test its existence.
if ~isfield(M_.(varexpectationmodel.auxiliary_model_type), varexpectationmodel.auxiliary_model_name)
    error('Unknown VAR (%s) in VAR_EXPECTATION_MODEL (%s)!', varexpectationmodel.auxiliary_model_name, varexpectationmodelname)
end

auxmodel = M_.(varexpectationmodel.auxiliary_model_type).(varexpectationmodel.auxiliary_model_name);

% Check that we have the values of the VAR matrices.
if ~isfield(oo_.(varexpectationmodel.auxiliary_model_type), varexpectationmodel.auxiliary_model_name)
    error('Auxiliary model %s has to be estimated or calibrated first!', varexpectationmodel.auxiliary_model_name)
end

auxcalib = oo_.(varexpectationmodel.auxiliary_model_type).(varexpectationmodel.auxiliary_model_name);

if ~isfield(auxcalib, 'CompanionMatrix') || any(isnan(auxcalib.CompanionMatrix(:)))
    message = sprintf('Auxiliary model %s has to be estimated first.', varexpectationmodel.auxiliary_model_name);
    message = sprintf('%s\nPlease use *get_companion_matrix command first.', message);
    error(message);
end

% Set discount factor
if isfield(varexpectationmodel, 'discount_value')
    discountfactor = varexpectationmodel.discount_value;
else
    if isfield(varexpectationmodel, 'discount_index')
        discountfactor = M_.params(varexpectationmodel.discount_index);
    else
        error('This is most likely a bug. Pleasse conntact the Dynare Team.')
    end
end

% A discount factor has to be positive.
if discountfactor<=0
    error('The discount factor must be positive.')
end

% A discount factor cannot be greater than one.
if discountfactor>1
    error('The discount cannot be greater than one.')
end

% Set timeshift
if isfield(varexpectationmodel, 'time_shift')
    timeshift = varexpectationmodel.time_shift;
else
    timeshift = 0;
end

% Throw an error if timeshift is positive
if timeshift>0
    error('Option time_shift of the var_expectation command cannot be positive.')
end

% Set variables_id in VAR model
m = length(varexpectationmodel.expr.vars);
variables_id_in_var = NaN(m,1);
for i = 1:m
    j = find(strcmp(auxmodel.list_of_variables_in_companion_var, M_.endo_names{varexpectationmodel.expr.vars(i)}));
    if isempty(j)
        error('Cannot find variable %s in the companion VAR', M_.endo_names{varexpectationmodel.expr.vars(i)})
    else
        variables_id_in_var(i) = find(strcmp(auxmodel.list_of_variables_in_companion_var, M_.endo_names{varexpectationmodel.expr.vars(i)}));
    end
end

if isfield(auxmodel, 'isconstant') && auxmodel.isconstant
    variables_id_in_var = variables_id_in_var+1;
end

% Get the horizon parameter.
horizon = varexpectationmodel.horizon;

% Check the horizon parameter
wrong_horizon_parameter = true;
if length(horizon)==1
    if isnumeric(horizon)
        if isfinite(horizon)
            if isint(horizon)
                if horizon>0
                    wrong_horizon_parameter = false;
                end
            end
        end
    end
elseif length(horizon)==2
    if isnumeric(horizon)
        if isfinite(horizon(1))
            if isint(horizon(1))
                if horizon(1)>=0
                    if isinf(horizon(2)) || (isint(horizon(2)) && horizon(2)>horizon(1)) 
                        wrong_horizon_parameter = false;
                    end
                end
            end
        end
    end
end

if wrong_horizon_parameter
    error('horizon must be an integer scalar or an integer vector with two elements.')
end

% Get the companion matrix
CompanionMatrix = auxcalib.CompanionMatrix;

% Get the dimension of the problem.
n = length(CompanionMatrix);

% Set the selection vector
alpha = zeros(1, n);
alpha(variables_id_in_var) = varexpectationmodel.expr.constants;
params_id_in_var = ~isnan(varexpectationmodel.expr.params);
alpha(variables_id_in_var(params_id_in_var)) = varexpectationmodel.expr.params(params_id_in_var);

if length(horizon)==1
    % Compute the reduced form parameters of the (discounted) forecast in period t+horizon(1)
    if varexpectationmodel.horizon==1
        parameters = discountfactor*(alpha*CompanionMatrix);
    elseif horizon>1 
        parameters = alpha*mpower(discountfactor*CompanionMatrix, varexpectationmodel.horizon);
    end
    if timeshift<0
        parameters = parameters*mpower(CompanionMatrix, -timeshift);
    end
else
    % Compute the reduced form parameters of the discounted sum of forecasts between t+horizon(1) and
    % t+horizon(2). Not that horzizon(2) need not be finite.
    if horizon(1)==0 && isinf(horizon(2))
        parameters = alpha/(eye(n)-discountfactor*CompanionMatrix);
        if timeshift<0
            parameters = parameters*mpower(CompanionMatrix, -timeshift);
        end
    elseif horizon(1)>0 && isinf(horizon(2))
        % Define the discounted companion matrix
        DiscountedCompanionMatrix = discountfactor*CompanionMatrix;
        % First compute the parameters implied by the discounted sum from h=0 to h=horizon(1)-1
        tmp1 = zeros(n);
        for h=0:horizon(1)-1
            tmp1 = tmp1 + mpower(DiscountedCompanionMatrix, h); 
        end
        tmp1 = alpha*tmp1;
        if timeshift<0
            tmp1 = tmp1*mpower(CompanionMatrix, -timeshift);
        end
        % Second compute the parameters implied by the discounted sum from h=0 to h=Inf 
        tmp2 = alpha/(eye(n)-DiscountedCompanionMatrix);
        if timeshift<0
            tmp2 = tmp2*mpower(CompanionMatrix, -timeshift);
        end
        % Finally
        parameters = tmp2-tmp1;
    elseif isfinite(horizon(2))
        % Define the discounted companion matrix
        DiscountedCompanionMatrix = discountfactor*CompanionMatrix;
        tmp = zeros(n);
        for h=horizon(1):horizon(2)
            tmp = tmp + mpower(DiscountedCompanionMatrix, h);
        end
        parameters = alpha*tmp;
        if timeshift<0
            parameters = parameters*mpower(CompanionMatrix, -timeshift);
        end
    end
end

% Update reduced form parameters in M_.params.
if isequal(varexpectationmodel.auxiliary_model_type, 'var')
    if M_.var.(varexpectationmodel.auxiliary_model_name).isconstant
        M_.params(varexpectationmodel.param_indices) = parameters;
    else
        M_.params(varexpectationmodel.param_indices(1)) = .0;
        M_.params(varexpectationmodel.param_indices(2:end)) = parameters;
    end
else
    M_.params(varexpectationmodel.param_indices) = parameters;
end