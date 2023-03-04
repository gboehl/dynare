function search(variablename, option)

% Prints equations where the variable appears in.
%
% INPUTS
% - variablename       [string]    Name of the variable to be traced.
% - option             [string]    if equal to 'withparamvalues' reports the values of the parameters, otherwise throws an error.
%
% OUTPUTS
% None

% Copyright © 2019-2023 Dynare Team
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

withexpansion = true;

if ~nargin
    error('Provide endogenous variable name as first input argument.');
end

jname = sprintf('%s%smodel%sjson%smodfile-original.json', M_.fname, filesep, filesep, filesep);

% Check if corresponding JSON file exists.
if exist(jname, 'file') ~= 2
    error('Please use the json option (See the Dynare invocation section in the reference manual).');
end

% Check that the first input is a character array.
if ~ischar(variablename)
    error('First input argument must be a string.');
end

% Check that the variable is actually a variable in the model.
if ~ismember(variablename, [M_.exo_names; M_.endo_names])
    error('There is no variable named %s!', variablename);
end

withparamvalues = false; % By default the values of the parameters are not printed

if nargin==2
    if ischar(option)
        if strcmp(option, 'withparamvalues')
            withparamvalues = true;
        else
            error('Unknown option %s.', option)
        end
    else
        error('Second input argument must be a row character array.')
    end
end

if nargin>2
    error('Too many input arguments.')
end

% Load the JSON file.
jsonfile = loadjson_(jname);
model = jsonfile.model;

% Print the equations the variable appears in.
for it = 1:length(M_.mapping.(variablename).eqidx)
    if M_.mapping.(variablename).eqidx(it)>length(model)
        % Equations appended by the preprocessor for auxiliary variables are not displayed, except if it is an aux variable
        % for the PAC expectation.
        % TODO Probably should do the same for VAR expectation.
        if isfield(M_, 'lhs') && length(M_.lhs{M_.mapping.(variablename).eqidx(it)})>15 && isequal(M_.lhs{M_.mapping.(variablename).eqidx(it)}(1:16), 'pac_expectation_')
            id = M_.mapping.(M_.lhs{M_.mapping.(variablename).eqidx(it)}).eqidx(1);
        else
            continue
        end
    else
        id = M_.mapping.(variablename).eqidx(it);
    end
    rhs = model{id}.rhs;
    if withexpansion
        if isfield(M_, 'pac') && contains(rhs, 'pac_expectation')
            % Get the index of the equation's PAC model.
            models = fieldnames(M_.pac);
            idx = find(~cellfun('isempty',cellfun(@(s)find(contains(rhs,s)),models,'uni',0)));
            % Get the expanded PAC_EXPECTATION term.
            [pac_expression, growthneutralitycorrection] = write_expectations(models{idx}, 'pac', true);
            expression = [sprintf('\n\t + %s', growthneutralitycorrection) TransformExpandedExpr(pac_expression)];
            rhs = strrep(rhs, ['+pac_expectation(model_name = ' models{idx} ')'], expression);
        elseif isfield(M_, 'var_expectation') && contains(rhs, 'var_expectation')
            % Get the index of the equation's VAR model.
            models = fieldnames(M_.var_expectation);
            idx = find(~cellfun('isempty',cellfun(@(s)find(contains(rhs,s)),models,'uni',0)));
            % Get the expanded VAR_EXPECTATION term.
            expression = write_expectations(models{idx}, 'var', true);
            expression = TransformExpandedExpr(expression);
            rhs = strrep(rhs, ['+var_expectation(model_name = ' models{idx} ')'], expression);
        elseif ~isfield(M_, 'pac') && ~isfield(M_, 'var_expectation')
            withexpansion = false;
        end
    end
    lhs = model{id}.lhs;
    if withparamvalues
        pnames = get_variables_and_parameters_in_equation(lhs, rhs, M_);
        if ~isempty(pnames)
            for ip = 1:length(pnames)
                iP = find(strcmp(M_.param_names, pnames{ip}));
                if isnan(M_.params(iP))
                    continue
                end
                value = M_.params(iP);
                rhs = regexprep(rhs, sprintf('\\<%s\\>', pnames{ip}), num2str(value, 6));
                lhs = regexprep(lhs, sprintf('\\<%s\\>', pnames{ip}), num2str(value, 6));
            end
        end
    end
    skipline()
    fprintf('%s = %s;\n', lhs, rhs);
end

skipline()

function [transformed_expression] = TransformExpandedExpr(expression)
transformed_expression = splitlines(expression);
transformed_expression{1} = sprintf(' + %s', transformed_expression{1});
transformed_expression = sprintf('\n\t%s', transformed_expression{:});
