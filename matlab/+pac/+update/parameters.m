function M_ = parameters(pacname, M_, oo_, verbose)
% M_ = parameters(pacname, M_, oo_, verbose)
% Updates the parameters of a PAC equation.
%
% INPUTS
% - pacname       [string]    Name of the pac equation.
% - M_            [struct]    M_ global structure (model properties)
% - oo_           [struct]    oo_ global structure (model results)
%
% OUTPUTS
% - none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2018-2024 Dynare Team
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

if nargin<4
    verbose = false;
end

% Check that the first input is a row character array.
if ~isrow(pacname)==1 || ~ischar(pacname)
    error('First input argument must be a row character array!')
end

% Check the name of the PAC model.
if ~isfield(M_.pac, pacname)
    error('PAC model %s is not defined in the model block!', pacname)
end

% Get PAC model description
pacmodel = M_.pac.(pacname);

if pacmodel.model_consistent_expectations
    error('This function cannot be used with Model Consistent Expectation. Try pac.mce.parameters instead.')
end

% Get the name of the associated auxiliary model (VAR or TREND_COMPONENT) model and test its existence.
if ~isfield(M_.(pacmodel.auxiliary_model_type), pacmodel.auxiliary_model_name)
    error('Unknown auxiliary model (%s) in PAC model (%s)!', pacmodel.auxiliary_model_name, pacname)
end
varmodel = M_.(pacmodel.auxiliary_model_type).(pacmodel.auxiliary_model_name);
% Check that we have the values of the VAR or TREND_COMPONENT matrices.
if ~isfield(oo_.(pacmodel.auxiliary_model_type), pacmodel.auxiliary_model_name)
    error('Auxiliary model %s has to be estimated first!', pacmodel.auxiliary_model_name)
end
varcalib = oo_.(pacmodel.auxiliary_model_type).(pacmodel.auxiliary_model_name);
if ~isfield(varcalib, 'CompanionMatrix') || any(isnan(varcalib.CompanionMatrix(:)))
    error('Auxiliary model %s has to be estimated first.', pacmodel.auxiliary_model_name)
end

% Show the equations where this PAC model is used.
if verbose
    fprintf('PAC model %s is used in equation %s.\n', pacname, pacmodel.eq_name);
    skipline()
end

% Do we need to decompose the PAC expectation?
if isfield(pacmodel, 'components')
    numberofcomponents = length(pacmodel.components);
else
    numberofcomponents = 0;
end

% Makes no sense to have a composite PAC target with a trend component auxiliary model (where all the variables are non stationnary)
if numberofcomponents>0 && strcmp(pacmodel.auxiliary_model_type, 'trend_component')
    error('Composite PAC target not allowed with trend component model.')
end

% Build the vector of PAC parameters (ECM parameter + autoregressive parameters).
pacvalues = M_.params([pacmodel.ec.params; pacmodel.ar.params(1:pacmodel.max_lag)']);

% Get the indices for the stationary/nonstationary variables in the VAR system.
if numberofcomponents
    id = cell(numberofcomponents, 1);
    for i=1:numberofcomponents
        id(i) = {find(strcmp(M_.endo_names{pacmodel.components(i).endo_var}, varmodel.list_of_variables_in_companion_var))};
        if isempty(id{i})
            % Find the auxiliary variables if any
            ad = find(cell2mat(cellfun(@(x) isauxiliary(x, [8 10]), varmodel.list_of_variables_in_companion_var, 'UniformOutput', false)));
            if isempty(ad)
                error('Cannot find the trend variable in the Companion VAR model.')
            else
                for j=1:length(ad)
                    auxinfo = M_.aux_vars(get_aux_variable_id(varmodel.list_of_variables_in_companion_var{ad(j)}));
                    if isequal(auxinfo.endo_index, pacmodel.components(i).endo_var)
                        id(i) = {ad(j)};
                        break
                    end
                    if isequal(auxinfo.type, 8) && isequal(auxinfo.orig_index, pacmodel.components(i).endo_var)
                        id(i) = {ad(j)};
                        break
                    end
                end
            end
            if isempty(id{i})
                error('Cannot find the trend variable in the Companion VAR model.')
            end
        end
    end
else
    id = {find(strcmp(M_.endo_names{pacmodel.ec.vars(pacmodel.ec.istarget)}, varmodel.list_of_variables_in_companion_var))};
    if isempty(id{1})
        % Find the auxiliary variables if any
        ad = find(cell2mat(cellfun(@(x) isauxiliary(x, [8 10]), varmodel.list_of_variables_in_companion_var, 'UniformOutput', false)));
        if isempty(ad)
            error('Cannot find the trend variable in the auxiliary VAR / Trend component model.')
        else
            for i=1:length(ad)
                auxinfo = M_.aux_vars(get_aux_variable_id(varmodel.list_of_variables_in_companion_var{ad(i)}));
                if isequal(auxinfo.endo_index, pacmodel.ec.vars(pacmodel.ec.istarget))
                    id = {ad(i)};
                    break
                end
                if isequal(auxinfo.type, 8) && isequal(auxinfo.orig_index, pacmodel.ec.vars(pacmodel.ec.istarget))
                    id = {ad(i)};
                    break
                end
            end
        end
        if isempty(id{1})
            error('Cannot find the trend variable in the auxiliary VAR / Trend component model.')
        end
    end
end

if ~numberofcomponents
    % Infer the kind of PAC exoectation
    if isequal(pacmodel.auxiliary_model_type, 'var')
        if varmodel.nonstationary(id{1})
            kind = {'dd'};
            if varmodel.isconstant
                id{1} = id{1}+1;
            end
        else
            kind = {'ll'};
            if varmodel.isconstant
                id{1} = id{1}+1;
            end
        end
    else
        % Trend component model is assumed.
        kind = {'dd'};
    end
else
    if varmodel.isconstant
        for i=1:numberofcomponents
            id{i} = id{i}+1;
        end
    end
end

% Override kind with the information provided by the user or update M_.pac
if ~numberofcomponents
    if ~isempty(pacmodel.kind)
        kind = {pacmodel.kind};
    else
        pacmodel.kind = kind{1};
    end
else
    kind = cell(numberofcomponents,1);
    for i=1:numberofcomponents
        if isempty(pacmodel.components(i).kind)
            error('kind declaration is mandatory for each component in pac_target_info.')
        else
            kind{i} = pacmodel.components(i).kind;
        end
    end
end

% Get the value of the discount factor.
beta = M_.params(pacmodel.discount_index);

% Is growth argument passed to pac_expectation?
if isfield(pacmodel, 'growth_str')
    growth_flag = true;
else
    growth_flag = false;
    for i=1:numberofcomponents
        if isfield(pacmodel.components(i), 'growth_str')
            growth_flag = true;
            break
        end
    end
end

% Do we have rule of thumb agents? γ is the share of optimizing agents.
if isfield(pacmodel, 'non_optimizing_behaviour')
    gamma = M_.params(pacmodel.share_of_optimizing_agents_index);
else
    gamma = 1.0;
end

% Get h vector (plus the parameter for the growth neutrality correction).
if growth_flag
    h = cell(1,length(id));
    growthneutrality = cell(1,length(id));
    for i=1:length(id)
        [h{i}, growthneutrality{i}] = hVectors([pacvalues; beta], varcalib.CompanionMatrix, pacmodel.auxiliary_model_type, kind{i}, id{i});
    end
else
    h = cell(1,length(id));
    for i=1:length(id)
        h(i) = {hVectors([pacvalues; beta], varcalib.CompanionMatrix, pacmodel.auxiliary_model_type, kind{i}, id{i})};
    end
end

% Update M_.params with h
if isequal(pacmodel.auxiliary_model_type, 'var')
    if M_.var.(pacmodel.auxiliary_model_name).isconstant
        if isfield(pacmodel, 'h_param_indices')
            % No decomposition
            M_.params(pacmodel.h_param_indices) = h{1};
        else
            for i=1:numberofcomponents
                M_.params(pacmodel.components(i).h_param_indices) = h{i};
            end
        end
    else
        if isfield(pacmodel, 'h_param_indices')
            % No decomposition
            M_.params(pacmodel.h_param_indices(1)) = .0;
            M_.params(pacmodel.h_param_indices(2:end)) = h{1};
        else
            for i=1:numberofcomponents
                M_.params(pacmodel.components(i).h_param_indices(1)) = .0;
                M_.params(pacmodel.components(i).h_param_indices(2:end)) = h{i};
            end
        end
    end % If the auxiliary model (VAR) has no constant.
else
    M_.params(pacmodel.h_param_indices) = h{1};
end % if auxiliary model is a VAR

% Update the parameter related to the growth neutrality correction.
if growth_flag
    % Growth neutrality as returned by hVector is valid iff
    % there is no exogenous variables in the model and in the
    % absence of non optimizing agents.
    for j=1:length(id)
        if isnan(growthneutrality{j})
            continue
        end
        gg = -(growthneutrality{j}-1); % Finite sum of autoregressive parameters + infinite sum of the coefficients in the PAC expectation term.
        cc = 1.0-gg*gamma;          % First adjustment of the growth neutrality correction (should also be divided by gamma, done below at the end of this section).
                                    % We may have to further change the correction if we have nonzero mean exogenous variables.
        ll = 0.0;
        if isfield(pacmodel, 'optim_additive')
            % Exogenous variables are present in the λ part (optimizing agents).
            tmp0 = 0;
            for i=1:length(pacmodel.optim_additive.params)
                if isnan(pacmodel.optim_additive.params(i)) && islogical(pacmodel.optim_additive.bgp{i}) && pacmodel.optim_additive.bgp{i}
                    tmp0 = tmp0 + pacmodel.optim_additive.scaling_factor(i);
                elseif ~isnan(pacmodel.optim_additive.params(i)) && islogical(pacmodel.optim_additive.bgp{i}) && pacmodel.optim_additive.bgp{i}
                    tmp0 = tmp0 + M_.params(pacmodel.optim_additive.params(i))*pacmodel.optim_additive.scaling_factor(i);
                elseif ~islogical(pacmodel.optim_additive.bgp{i})
                    error('It is not possible to provide a value for the mean of an exogenous variable appearing in the optimal part of the PAC equation.')
                end
            end
            cc = cc - tmp0*gamma;
        end
        if gamma<1
            if isfield(pacmodel, 'non_optimizing_behaviour') && isfield(pacmodel.non_optimizing_behaviour, 'params')
                % Exogenous variables are present in the 1-λ part (rule of thumb agents).
                tmp0 = 0;
                tmp1 = 0;
                for i=1:length(pacmodel.non_optimizing_behaviour.params)
                    if isnan(pacmodel.non_optimizing_behaviour.params(i)) && islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && pacmodel.non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + pacmodel.non_optimizing_behaviour.scaling_factor(i);
                    elseif ~isnan(pacmodel.non_optimizing_behaviour.params(i)) && islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && pacmodel.non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + M_.params(pacmodel.non_optimizing_behaviour.params(i))*pacmodel.non_optimizing_behaviour.scaling_factor(i);
                    elseif ~islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && isnumeric(pacmodel.non_optimizing_behaviour.bgp{i}) && isnan(pacmodel.non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + pacmodel.non_optimizing_behaviour.scaling_factor(i)*pacmodel.non_optimizing_behaviour.bgp{i};
                    elseif ~islogical(pacmodel.non_optimizing_behaviour.bgp{i}) && isnumeric(pacmodel.non_optimizing_behaviour.bgp{i}) && ~isnan(pacmodel.non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + pacmodel.non_optimizing_behaviour.scaling_factor(i)*pacmodel.non_optimizing_behaviour.params(i)*pacmodel.non_optimizing_behaviour.bgp{i};
                    end
                end
                cc = cc - (1.0-gamma)*tmp0;
                ll = -(1.0-gamma)*tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
            end
        end
        if isfield(pacmodel, 'additive')
            % Exogenous variables are present outside of the λ and (1-λ) parts (or we have exogenous variables in a "pure" PAC equation.
            tmp0 = 0;
            tmp1 = 0;
            for i=1:length(pacmodel.additive.params)
                if isnan(pacmodel.additive.params(i)) && islogical(pacmodel.additive.bgp{i}) && pacmodel.additive.bgp{i}
                    tmp0 = tmp0 + pacmodel.additive.scaling_factor(i);
                elseif ~isnan(pacmodel.additive.params(i)) && islogical(pacmodel.additive.bgp{i}) && pacmodel.additive.bgp{i}
                    tmp0 = tmp0 + M_.params(pacmodel.additive.params(i))*pacmodel.additive.scaling_factor(i);
                elseif ~islogical(pacmodel.additive.bgp{i}) && isnumeric(pacmodel.additive.bgp{i}) && isnan(pacmodel.additive.params(i))
                    tmp1 = tmp1 + pacmodel.additive.scaling_factor(i)*pacmodel.additive.bgp{i};
                elseif ~islogical(pacmodel.additive.bgp{i}) && isnumeric(pacmodel.additive.bgp{i}) && ~isnan(pacmodel.additive.params(i))
                    tmp1 = tmp1 + pacmodel.additive.scaling_factor(i)*pacmodel.additive.params(i)*pacmodel.additive.bgp{i};
                end
            end
            cc = cc - tmp0;
            ll = ll - tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
        end
        if isfield(pacmodel, 'growth_neutrality_param_index')
            M_.params(pacmodel.growth_neutrality_param_index) = cc/gamma; % Multiplies the variable or expression provided though the growth option in command pac_model.
        else
            M_.params(pacmodel.components(j).growth_neutrality_param_index) = cc/gamma; % Multiplies the variable or expression provided though the growth option in command pac_model.
        end
    end
end
