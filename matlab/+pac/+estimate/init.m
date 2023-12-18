function [pacmodl, lhs, rhs, pnames, enames, xnames, rname, pid, eid, xid, pnames_, ipnames_, params, data, islaggedvariables] =  init(M_, oo_, eqname, params, data, range)

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

% Get the original equation to be estimated
[~, RHS] = get_lhs_and_rhs(eqname, M_, true);

% Check that the equation has a PAC expectation term.
if ~contains(RHS, 'pac_expectation', 'IgnoreCase', true)
    error('This is not a PAC equation.')
end

% Get the name of the PAC model.
pattern = '(\(model_name\s*=\s*)(?<name>\w+)\)';
pacmodl = regexp(RHS, pattern, 'names');
pacmodl = pacmodl.name;

pacmodel = M_.pac.(pacmodl);

% Get the transformed equation to be estimated.
[lhs, rhs, json] = get_lhs_and_rhs(eqname, M_);

% Get definition of aux. variable for pac expectation...
if isfield(pacmodel, 'aux_id')
    auxrhs = {M_.lhs{pacmodel.aux_id}, json.model{pacmodel.aux_id}.rhs};
elseif isfield(pacmodel, 'components')
    auxrhs = cell(length(pacmodel.components), 2);
    for i=1:length(pacmodel.components)
        auxrhs{i,1} = M_.lhs{pacmodel.components(i).aux_id};
        auxrhs{i,2} = sprintf('(%s)', json.model{pacmodel.components(i).aux_id}.rhs);
    end
else
    error('Cannot find auxiliary variables for PAC expectation.')
end

% ... and substitute in rhs.
for i=1:rows(auxrhs)
    rhs = strrep(rhs, auxrhs{i,1}, auxrhs{i,2});
end

% Get pacmodel properties
pacmodel = M_.pac.(pacmodl);

% Get the parameters and variables in the PAC equation.
[pnames, enames, xnames, pid, eid, xid] = get_variables_and_parameters_in_equation(lhs, rhs, M_);

% Get list and indices of estimated parameters.
ipnames_ = get_estimated_parameters_indices(params, pnames, eqname, M_);

% If equation is estimated by recursive OLS, ensure that the error
% correction parameter comes first, followed by the autoregressive
% parameters (in increasing order w.r.t. the lags).
stack = dbstack;
ipnames__ = ipnames_;                              % The user provided order.
if isequal(stack(2).name, 'iterative_ols')
    ipnames_ = [pacmodel.ec.params; pacmodel.ar.params'];
    if isfield(pacmodel, 'optim_additive')
        ipnames_ = [ipnames_; pacmodel.optim_additive.params(~isnan(pacmodel.optim_additive.params))'];
    end
    if isfield(pacmodel, 'additive')
        ipnames_ = [ipnames_; pacmodel.additive.params(~isnan(pacmodel.additive.params))'];
    end
    if isfield(pacmodel, 'share_of_optimizing_agents_index')
        ipnames_ = [ipnames_; pacmodel.share_of_optimizing_agents_index];
    end
    for i=1:length(ipnames_)
        if ~ismember(ipnames_(i), ipnames__)
            % This parameter is not estimated.
            ipnames_(i) = NaN;
        end
    end
end

% Remove calibrated parameters (if any).
ipnames_(isnan(ipnames_)) = [];

% Reorder params if needed.
[~, permutation] = ismember(ipnames__, ipnames_);
pnames_ = fieldnames(params);
pnames_ = pnames_(permutation);
params = orderfields(params, permutation);

% Add the auxiliary variables in the dataset.
if M_.endo_nbr>M_.orig_endo_nbr
    data = feval([M_.fname '.dynamic_set_auxiliary_series'], data, M_.params);
end

% Check that the data for endogenous variables have values.
if any(isnan(data{enames{:}}(range).data(:)))
    error('Some variable values are missing in the database.')
end

% Set the number of exogenous variables.
xnbr = length(xnames);

% Test if we have a residual and get its name (-> rname).
if isequal(xnbr, 1)
    rname = M_.exo_names{strcmp(xnames{1}, M_.exo_names)};
    if ~all(isnan(data{xnames{1}}.data))
        error('The residual (%s) must have NaN values in the provided database.', xnames{1})
    end
else
    % We have observed exogenous variables in the PAC equation.
    tmp = data{xnames{:}}(range).data;
    idx = find(all(~isnan(tmp))); % Indices for the observed exogenous variables.
    if isequal(length(idx), length(xnames))
        error('There is no residual in this equation, all the exogenous variables are observed.')
    else
        if length(idx)<length(xnames)-1
            error('It is not allowed to have more than one residual in a PAC equation')
        end
        irname = setdiff(1:length(xnames), idx);
        rname = xnames{irname};
    end
end


% Remove residuals from the equation.
%
% Note that a plus or minus will remain in the equation, but this seems to
% be without consequence.
rhs = regexprep(rhs, rname, '');

% Create a dummy variable
islaggedvariables = contains(rhs, '(-1)');