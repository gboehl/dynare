function [expression, growthneutralitycorrection] = write_expectations(eqname, expectationmodelname, expectationmodelkind, iscrlf)

% Prints the exansion of the VAR_EXPECTATION or PAC_EXPECTATION term in files.
%
% INPUTS
% - eqname                      [string]    Name of the equation.
% - epxpectationmodelname       [string]    Name of the expectation model.
% - expectationmodelkind        [string]    Kind of the expectation model ('var' or 'pac').
% - iscrlf                      [string]    Adds carriage return after each additive term if true.
%
% OUTPUTS
% - expression                  [string]    Unrolled expectation expression.
% - growthneutralitycorrection  [string]

% Copyright © 2019-2021 Dynare Team
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

global M_

if ismember(expectationmodelkind, {'var', 'pac'})
    if isequal(expectationmodelkind, 'var')
        expectationmodelfield = 'var_expectation';
    else
        expectationmodelfield = 'pac';
        % Get the equation tag (in M_.pac.(pacmodl).equations)
        eqtag = M_.pac.(expectationmodelname).tag_map{strcmp(M_.pac.(expectationmodelname).tag_map(:,1), eqname),2};
    end
else
    error('Value of third input argument must be ''var'' or ''pac''.')
end

expectationmodel = M_.(expectationmodelfield).(expectationmodelname);

if nargout>1 && isequal(expectationmodelkind, 'var')
    error('Cannot return more than one argument if the expectation model is a VAR.')
end

if nargin<4
    iscrlf = false;
end

% Get the name of the associated VAR model and test its existence.
if ~isfield(M_.(expectationmodel.auxiliary_model_type), expectationmodel.auxiliary_model_name)
    switch expectationmodelkind
      case 'var-expectations'
        error('Unknown VAR/TREND_COMPONENT model (%s) in VAR_EXPECTATION_MODEL (%s)!', expectationmodel.auxiliary_model_name, expectationmodelname)
      case 'pac-expectations'
        error('Unknown VAR/TREND_COMPONENT model (%s) in PAC_EXPECTATION_MODEL (%s)!', expectationmodel.auxiliary_model_name, expectationmodelname)
      otherwise
    end
end

auxmodel = M_.(expectationmodel.auxiliary_model_type).(expectationmodel.auxiliary_model_name);

maxlag = max(auxmodel.max_lag);
if isequal(expectationmodel.auxiliary_model_type, 'trend_component')
    % Need to add a lag since the error correction equations are rewritten in levels.
    maxlag = maxlag+1;
end

id = 0;

if isequal(expectationmodelkind, 'var')
    timeindices = (0:(maxlag-1))+abs(expectationmodel.time_shift);
end

if isequal(expectationmodelkind, 'var') && isequal(expectationmodel.auxiliary_model_type, 'var')
    id = id+1;
    expression = sprintf('%s', M_.param_names{expectationmodel.param_indices(id)});
end

if isequal(expectationmodelkind, 'pac') && isequal(expectationmodel.auxiliary_model_type, 'var')
    id = id+1;
    expression = sprintf('%s+%s', M_.param_names{expectationmodel.equations.(eqtag).h0_param_indices(id)}, ...
                         M_.param_names{expectationmodel.equations.(eqtag).h1_param_indices(id)});
end

for i=1:maxlag
    for j=1:length(auxmodel.list_of_variables_in_companion_var)
        id = id+1;
        variable = auxmodel.list_of_variables_in_companion_var{j};
        transformations = {};
        ida = get_aux_variable_id(variable);
        op = 0;
        while ida
            op = op+1;
            if isequal(M_.aux_vars(ida).type, 8)
                transformations(op) = {'diff'};
                variable = M_.endo_names{M_.aux_vars(ida).orig_index};
                ida = get_aux_variable_id(variable);
            elseif isequal(M_.aux_vars(ida).type, 10)
                transformations(op) = {M_.aux_vars(ida).unary_op};
                variable = M_.endo_names{M_.aux_vars(ida).orig_index};
                ida = get_aux_variable_id(variable);
            else
                error('This case is not implemented.')
            end
        end
        switch expectationmodelkind
          case 'var'
            parameter = M_.param_names{expectationmodel.param_indices(id)};
          case 'pac'
            parameter = '';
            if ~isempty(expectationmodel.equations.(eqtag).h0_param_indices)
                parameter = M_.param_names{expectationmodel.equations.(eqtag).h0_param_indices(id)};
            end
            if ~isempty(expectationmodel.equations.(eqtag).h1_param_indices)
                if isempty(parameter)
                    parameter = M_.param_names{expectationmodel.equations.(eqtag).h1_param_indices(id)};
                else
                    parameter = sprintf('(%s+%s)', parameter, M_.param_names{expectationmodel.equations.(eqtag).h1_param_indices(id)});
                end
            end
          otherwise
        end
        switch expectationmodelkind
          case 'var'
            if timeindices(i)
                variable = sprintf('%s(-%d)', variable, timeindices(i));
            end
          case 'pac'
            variable = sprintf('%s(-%d)', variable, i);
          otherwise
        end
        if ~isempty(transformations)
            for k=length(transformations):-1:1
                variable = sprintf('%s(%s)', transformations{k}, variable);
            end
        end
        if isequal(id, 1)
            if iscrlf
                expression = sprintf('%s*%s\n', parameter, variable);
            else
                expression = sprintf('%s*%s', parameter, variable);
            end
        else
            if iscrlf
                expression = sprintf('%s + %s*%s\n', expression, parameter, variable);
            else
                expression = sprintf('%s + %s*%s', expression, parameter, variable);
            end
        end
    end
end

if isfield(expectationmodel, 'growth_neutrality_param_index')
    if numel(expectationmodel.growth_linear_comb) == 1
        growthneutralitycorrection = sprintf('%s*%s', M_.param_names{expectationmodel.growth_neutrality_param_index}, expectationmodel.growth_str);
    else
        growthneutralitycorrection = sprintf('%s*(%s)', M_.param_names{expectationmodel.growth_neutrality_param_index}, expectationmodel.growth_str);
    end
else
    growthneutralitycorrection = '';
end

if nargout==1 && ~isempty(growthneutralitycorrection)
    expression = sprintf('%s + %s', expression, growthneutralitycorrection);
end