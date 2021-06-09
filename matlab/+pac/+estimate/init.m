function [pacmodl, lhs, rhs, pnames, enames, xnames, rname, pid, eid, xid, pnames_, ipnames_, params, data, islaggedvariables, eqtag] = ...
    init(M_, oo_, eqname, params, data, range)

% Copyright (C) 2018-2019 Dynare Team
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
[LHS, RHS] = get_lhs_and_rhs(eqname, M_, true);

% Check that the equation has a PAC expectation term.
if ~contains(RHS, 'pac_expectation', 'IgnoreCase', true)
    error('This is not a PAC equation.')
end

% Get the name of the PAC model.
pattern = '(\(model_name\s*=\s*)(?<name>\w+)\)';
pacmodl = regexp(RHS, pattern, 'names');
pacmodl = pacmodl.name;

% Get the transformed equation to be estimated.
[lhs, rhs] = get_lhs_and_rhs(eqname, M_);

% Get the equation tag (in M_.pac.(pacmodl).equations)
eqtag = M_.pac.(pacmodl).tag_map{strcmp(M_.pac.(pacmodl).tag_map(:,1), eqname),2};

% Get the parameters and variables in the PAC equation.
[pnames, enames, xnames, pid, eid, xid] = get_variables_and_parameters_in_equation(lhs, rhs, M_);

% Get the list of estimated parameters
pnames_ = fieldnames(params);

% Check that the estimated parameters are used in the PAC equation.
ParametersNotInPAC = setdiff(pnames_, pnames);
if ~isempty(ParametersNotInPAC)
    skipline()
    if length(ParametersNotInPAC)>1
        list = sprintf('  %s\n', ParametersNotInPAC{:});
        remk = sprintf('  The following parameters:\n\n%s\n  do not appear in the PAC equation.', list);
    else
        remk = sprintf('  Parameter %s does not appear in the PAC equation.', ParametersNotInPAC{1});
    end
    disp(remk)
    skipline()
    error('The estimated parameters must be used in equation %s.', eqname)
end

% Get indices of estimated parameters.
ipnames_ = zeros(size(pnames_));
for i=1:length(ipnames_)
    ipnames_(i) = find(strcmp(pnames_{i}, M_.param_names));
end

% If equation is estimated by recursive OLS, ensure that the error
% correction parameter comes first, followed by the autoregressive
% parameters (in increasing order w.r.t. the lags).
stack = dbstack;
ipnames__ = ipnames_;                              % The user provided order.
if isequal(stack(2).name, 'iterative_ols')
    ipnames_ = [M_.pac.(pacmodl).equations.(eqtag).ec.params; M_.pac.(pacmodl).equations.(eqtag).ar.params'];
    if isfield(M_.pac.(pacmodl).equations.(eqtag), 'optim_additive')
        ipnames_ = [ipnames_; M_.pac.(pacmodl).equations.(eqtag).optim_additive.params(~isnan(M_.pac.(pacmodl).equations.(eqtag).optim_additive.params))'];
    end
    if isfield(M_.pac.(pacmodl).equations.(eqtag), 'additive')
        ipnames_ = [ipnames_; M_.pac.(pacmodl).equations.(eqtag).additive.params(~isnan(M_.pac.(pacmodl).equations.(eqtag).additive.params))'];
    end
    if isfield(M_.pac.(pacmodl).equations.(eqtag), 'share_of_optimizing_agents_index')
        ipnames_ = [ipnames_; M_.pac.(pacmodl).equations.(eqtag).share_of_optimizing_agents_index];
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
pnames_ = pnames_(permutation);
params = orderfields(params, permutation);

% Add the auxiliary variables in the dataset.
data = feval([M_.fname '.dynamic_set_auxiliary_series'], data, M_.params);

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