function DynareModel = parameters(pacname, DynareModel, DynareOutput, verbose)

% Updates the parameters of a PAC equation.
%
% INPUTS
% - pacname       [string]    Name of the pac equation.
% - DynareModel   [struct]    M_ global structure (model properties)
% - DynareOutput  [struct]    oo_ global structure (model results)
%
% OUTPUTS
% - none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2018-2021 Dynare Team
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
    verbose = true;
end

% Check that the first input is a row character array.
if ~isrow(pacname)==1 || ~ischar(pacname)
    error('First input argument must be a row character array!')
end

% Check the name of the PAC model.
if ~isfield(DynareModel.pac, pacname)
    error('PAC model %s is not defined in the model block!', pacname)
end

% Get PAC model description
pacmodel = DynareModel.pac.(pacname);

if pacmodel.model_consistent_expectations
    error('This function cannot be used with Model Consistent Expectation. Try pac.mce.parameters instead.')
end

% Get the name of the associated auxiliary model (VAR or TREND_COMPONENT) model and test its existence.
if ~isfield(DynareModel.(pacmodel.auxiliary_model_type), pacmodel.auxiliary_model_name)
    error('Unknown auxiliary model (%s) in PAC model (%s)!', pacmodel.auxiliary_model_name, pacname)
end
varmodel = DynareModel.(pacmodel.auxiliary_model_type).(pacmodel.auxiliary_model_name);
% Check that we have the values of the VAR or TREND_COMPONENT matrices.
if ~isfield(DynareOutput.(pacmodel.auxiliary_model_type), pacmodel.auxiliary_model_name)
    error('Auxiliary model %s has to be estimated first!', pacmodel.auxiliary_model_name)
end
varcalib = DynareOutput.(pacmodel.auxiliary_model_type).(pacmodel.auxiliary_model_name);
if ~isfield(varcalib, 'CompanionMatrix') || any(isnan(varcalib.CompanionMatrix(:)))
    error('Auxiliary model %s has to be estimated first.', pacmodel.auxiliary_model_name)
end

% Show the equations where this PAC model is used.
number_of_pac_eq = size(pacmodel.tag_map, 1);
if verbose
    fprintf('PAC model %s is used in %u equation(s):\n', pacname, number_of_pac_eq);
    skipline()
    for i=1:number_of_pac_eq
        fprintf('    - %s\n', pacmodel.tag_map{i,1});
    end
    skipline()
end

equations = pacmodel.equations;

for e=1:number_of_pac_eq
    eqtag = pacmodel.tag_map{e,2};
    % Build the vector of PAC parameters (ECM parameter + autoregressive parameters).
    pacvalues = DynareModel.params([equations.(eqtag).ec.params; equations.(eqtag).ar.params(1:equations.(eqtag).max_lag)']);
    % Get the indices for the stationary/nonstationary variables in the VAR system.
    id = find(strcmp(DynareModel.endo_names{equations.(eqtag).ec.vars(equations.(eqtag).ec.istarget)}, varmodel.list_of_variables_in_companion_var));
    if isempty(id)
        % Find the auxiliary variables if any
        ad = find(cell2mat(cellfun(@(x) isauxiliary(x, [8 10]), varmodel.list_of_variables_in_companion_var, 'UniformOutput', false)));
        if isempty(ad)
            error('Cannot find the trend variable in the Companion VAR/VECM model.')
        else
            for i=1:length(ad)
                auxinfo = DynareModel.aux_vars(get_aux_variable_id(varmodel.list_of_variables_in_companion_var{ad(i)}));
                if isequal(auxinfo.endo_index, equations.(eqtag).ec.vars(equations.(eqtag).ec.istarget))
                    id = ad(i);
                    break
                end
                if isequal(auxinfo.type, 8) && isequal(auxinfo.orig_index, equations.(eqtag).ec.vars(equations.(eqtag).ec.istarget))
                    id = ad(i);
                    break
                end
            end
        end
        if isempty(id)
            error('Cannot find the trend variable in the Companion VAR/VECM model.')
        end
    end
    if isequal(pacmodel.auxiliary_model_type, 'var')
        if varmodel.nonstationary(id)
            idns = id;
            ids = [];
            if varmodel.isconstant
                idns = idns+1;
            end
        else
            idns = [];
            ids = id;
            if varmodel.isconstant
                ids = ids+1;
            end
        end
    else
        % Trend component model is assumed.
        ids = [];
        idns = id;
    end
    % Get the value of the discount factor.
    beta = DynareModel.params(pacmodel.discount_index);
    % Is growth argument passed to pac_expectation?
    if isfield(pacmodel, 'growth_str')
        growth_flag = true;
    else
        growth_flag = false;
    end
    % Do we have rule of thumb agents? γ is the share of optimizing agents.
    if isfield(equations.(eqtag), 'non_optimizing_behaviour')
        gamma = DynareModel.params(equations.(eqtag).share_of_optimizing_agents_index);
    else
        gamma = 1.0;
    end
    % Get h0 and h1 vectors (plus the parameter for the growth neutrality correction).
    if growth_flag
        [h0, h1, growthneutrality] = hVectors([pacvalues; beta], varcalib.CompanionMatrix, ids, idns, pacmodel.auxiliary_model_type);
    else
        [h0, h1] = hVectors([pacvalues; beta], varcalib.CompanionMatrix, ids, idns, pacmodel.auxiliary_model_type);
    end
    % Update the parameters related to the stationary components.
    if ~isempty(h0)
        if isequal(pacmodel.auxiliary_model_type, 'var')
            if DynareModel.var.(pacmodel.auxiliary_model_name).isconstant
                DynareModel.params(equations.(eqtag).h0_param_indices) = h0;
            else
                DynareModel.params(equations.(eqtag).h0_param_indices(1)) = .0;
                DynareModel.params(equations.(eqtag).h0_param_indices(2:end)) = h0;
            end
        else
            DynareModel.params(equations.(eqtag).h0_param_indices) = h0;
        end
        DynareModel.params(pacmodel.equations.(eqtag).h0_param_indices) = h0;
    else
        if ~isempty(equations.(eqtag).h0_param_indices)
            DynareModel.params(equations.(eqtag).h0_param_indices) = .0;
        end
    end
    % Update the parameters related to the nonstationary components.
    if ~isempty(h1)
        if isequal(pacmodel.auxiliary_model_type, 'var')
            if DynareModel.var.(pacmodel.auxiliary_model_name).isconstant
                DynareModel.params(equations.(eqtag).h1_param_indices) = h1;
            else
                DynareModel.params(equations.(eqtag).h1_param_indices(1)) = .0;
                DynareModel.params(equations.(eqtag).h1_param_indices(2:end)) = h1;
            end
        else
            DynareModel.params(equations.(eqtag).h1_param_indices) = h1;
        end
    else
        if ~isempty(equations.(eqtag).h1_param_indices)
            DynareModel.params(equations.(eqtag).h1_param_indices) = .0;
        end
    end
    % Update the parameter related to the growth neutrality correction.
    if growth_flag
        % Growth neutrality as returned by hVector is valid iff
        % there is no exogenous variables in the model and in the
        % absence of non optimizing agents.
        gg = -(growthneutrality-1); % Finite sum of autoregressive parameters + infinite sum of the coefficients in the PAC expectation term.
        cc = 1.0-gg*gamma;          % First adjustment of the growth neutrality correction (should also be divided by gamma, done below at the end of this section).
                                    % We may have to further change the correction if we have nonzero mean exogenous variables.
        ll = 0.0;
        if isfield(equations.(eqtag), 'optim_additive')
            % Exogenous variables are present in the λ part (optimizing agents).
            tmp0 = 0;
            for i=1:length(equations.(eqtag).optim_additive.params)
                if isnan(equations.(eqtag).optim_additive.params(i)) && islogical(equations.(eqtag).optim_additive.bgp{i}) && equations.(eqtag).optim_additive.bgp{i}
                    tmp0 = tmp0 + equations.(eqtag).optim_additive.scaling_factor(i);
                elseif ~isnan(equations.(eqtag).optim_additive.params(i)) && islogical(equations.(eqtag).optim_additive.bgp{i}) && equations.(eqtag).optim_additive.bgp{i}
                    tmp0 = tmp0 + DynareModel.params(equations.(eqtag).optim_additive.params(i))*equations.(eqtag).optim_additive.scaling_factor(i);
                elseif ~islogical(equations.(eqtag).optim_additive.bgp{i})
                    error('It is not possible to provide a value for the mean of an exogenous variable appearing in the optimal part of the PAC equation.')
                end
            end
            cc = cc - tmp0*gamma;
        end
        if gamma<1
            if isfield(equations.(eqtag), 'non_optimizing_behaviour') && isfield(equations.(eqtag).non_optimizing_behaviour, 'params')
                % Exogenous variables are present in the 1-λ part (rule of thumb agents).
                tmp0 = 0;
                tmp1 = 0;
                for i=1:length(equations.(eqtag).non_optimizing_behaviour.params)
                    if isnan(equations.(eqtag).non_optimizing_behaviour.params(i)) && islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && equations.(eqtag).non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i);
                    elseif ~isnan(equations.(eqtag).non_optimizing_behaviour.params(i)) && islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && equations.(eqtag).non_optimizing_behaviour.bgp{i}
                        tmp0 = tmp0 + DynareModel.params(equations.(eqtag).non_optimizing_behaviour.params(i))*equations.(eqtag).non_optimizing_behaviour.scaling_factor(i);
                    elseif ~islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnumeric(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnan(equations.(eqtag).non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i)*equations.(eqtag).non_optimizing_behaviour.bgp{i};
                    elseif ~islogical(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && isnumeric(equations.(eqtag).non_optimizing_behaviour.bgp{i}) && ~isnan(equations.(eqtag).non_optimizing_behaviour.params(i))
                        tmp1 = tmp1 + equations.(eqtag).non_optimizing_behaviour.scaling_factor(i)*equations.(eqtag).non_optimizing_behaviour.params(i)*equations.(eqtag).non_optimizing_behaviour.bgp{i};
                    end
                end
                cc = cc - (1.0-gamma)*tmp0;
                ll = -(1.0-gamma)*tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
            end
        end
        if isfield(equations.(eqtag), 'additive')
            % Exogenous variables are present outside of the λ and (1-λ) parts (or we have exogenous variables in a "pure" PAC equation.
            tmp0 = 0;
            tmp1 = 0;
            for i=1:length(equations.(eqtag).additive.params)
                if isnan(equations.(eqtag).additive.params(i)) && islogical(equations.(eqtag).additive.bgp{i}) && equations.(eqtag).additive.bgp{i}
                    tmp0 = tmp0 + equations.(eqtag).additive.scaling_factor(i);
                elseif ~isnan(equations.(eqtag).additive.params(i)) && islogical(equations.(eqtag).additive.bgp{i}) && equations.(eqtag).additive.bgp{i}
                    tmp0 = tmp0 + DynareModel.params(equations.(eqtag).additive.params(i))*equations.(eqtag).additive.scaling_factor(i);
                elseif ~islogical(equations.(eqtag).additive.bgp{i}) && isnumeric(equations.(eqtag).additive.bgp{i}) && isnan(equations.(eqtag).additive.params(i))
                    tmp1 = tmp1 + equations.(eqtag).additive.scaling_factor(i)*equations.(eqtag).additive.bgp{i};
                elseif ~islogical(equations.(eqtag).additive.bgp{i}) && isnumeric(equations.(eqtag).additive.bgp{i}) && ~isnan(equations.(eqtag).additive.params(i))
                    tmp1 = tmp1 + equations.(eqtag).additive.scaling_factor(i)*equations.(eqtag).additive.params(i)*equations.(eqtag).additive.bgp{i};
                end
            end
            cc = cc - tmp0;
            ll = ll - tmp1/gamma; % TODO: ll should be added as a constant in the PAC equation (under the λ part) when unrolling pac_expectation.
        end
        DynareModel.params(pacmodel.growth_neutrality_param_index) = cc/gamma; % Multiplies the variable or expression provided though the growth option in command pac_model.
    end
end