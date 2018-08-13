function DynareModel = parameters(pacname, DynareModel, DynareOutput)

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

% Get the name of the associated VAR model and test its existence.
if ~isfield(DynareModel.var, pacmodel.var_model_name)
    error('Unknown VAR (%s) in PAC model (%s)!', pacmodel.var_model_name, pacname)
end

varmodel = DynareModel.var.(pacmodel.var_model_name);

% Check that we have the values of the VAR matrices.
if ~isfield(DynareOutput.var, pacmodel.var_model_name) 
    error('VAR model %s has to be estimated first!', pacmodel.var_model_name)
end

varcalib = DynareOutput.var.(pacmodel.var_model_name);

if ~isfield(varcalib, 'CompanionMatrix') || any(isnan(varcalib.CompanionMatrix(:)))
    error('VAR model %s has to be estimated first.', pacmodel.var_model_name)
end

% Build the vector of PAC parameters (ECM parameter + autoregressive parameters).
pacvalues = DynareModel.params([pacmodel.ec.params; pacmodel.ar.params(1:pacmodel.max_lag)']);

% Get the indices for the stationary/nonstationary variables in the VAR system.
id = find(strcmp(DynareModel.endo_names{pacmodel.ec.vars(1)}, varmodel.list_of_variables_in_companion_var));

if isempty(id)
    % Find the auxiliary variables if any
    ad = find(cell2mat(cellfun(@(x) isauxiliary(x, 8), varmodel.list_of_variables_in_companion_var, 'UniformOutput', false)));
    if isempty(ad)
        error('Cannot find the trend variable in the Companion VAR/VECM model.')
    else
        for i=1:length(ad)
            auxinfo = DynareModel.aux_vars(get_aux_variable_id(varmodel.list_of_variables_in_companion_var{ad(i)}));
            if isequal(auxinfo.orig_index, pacmodel.ec.vars(1))
                id = ad(i);
                break
            end
        end
    end
    if isempty(id)
        error('Cannot find the trend variable in the Companion VAR/VECM model.')
    end
end

if varmodel.nonstationary(id)
    idns = id;
    ids = [];
else
    idns = [];
    ids = id;
end

% Get the value of the discount factor.
beta = DynareModel.params(pacmodel.discount_index);

% Is growth argument passed to pac_expectation?
if isfield(pacmodel, 'growth_index')
    growth_flag = true;
    growth_type = pacmodel.growth_type;
else
    growth_flag = false;
end

% Get h0 and h1 vectors (plus the parameter for the growth neutrality correction).
if growth_flag
    [h0, h1, growthneutrality] = hVectors([pacvalues; beta], varcalib.CompanionMatrix, ids, idns, pacmodel.auxmodel);
else
    [h0, h1] = hVectors([pacvalues; beta], varcalib.CompanionMatrix, ids, idns, pacmodel.auxmodel);
end

% Update the parameters related to the stationary components.
if length(h0)
    DynareModel.params(pacmodel.h0_param_indices) = h0;
else
    if isfield(pacmodel, 'h0_param_indices') && length(pacmodel.h0_param_indices)
        DynareModel.params(pacmodel.h0_param_indices) = .0;
    end
end

% Update the parameters related to the nonstationary components.
if length(h1)
    DynareModel.params(pacmodel.h1_param_indices) = h1;
else
    if isfield(pacmodel, 'h1_param_indices') && length(pacmodel.h1_param_indices)
        DynareModel.params(pacmodel.h1_param_indices) = .0;
    end
end

% Update the parameter related to the growth neutrality correction.
if growth_flag
    DynareModel.params(pacmodel.growth_neutrality_param_index) = growthneutrality;
end