function nls(eqname, params, data, range)

% Estimates the parameters of a PAC equation by Nonlinear Least Squares.
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

% Get the parameters and variables in the PAC equation.
[pnames, enames, xnames, pid, eid, xid] = get_variables_and_parameters_in_equation(lhs, rhs, M_);

% List of objects to be replaced
objNames = [pnames; enames; xnames];
objIndex = [pid; eid; xid];
objTypes = [ones(length(pid), 1); 2*ones(length(eid), 1); 3*ones(length(eid), 1);];

[~,I] = sort(cellfun(@length, objNames), 'descend');
objNames = objNames(I);
objIndex = objIndex(I);
objTypes = objTypes(I);

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
    ipnames_(i) = strmatch(pnames_{i}, M_.param_names, 'exact');
end

% Add the auxiliary variables in the dataset.
data = feval([M_.fname '.dynamic_set_auxiliary_series'], data, M_.params);

% Check that the data for endogenous variables have values. Note that we
% also need data in range(1)-1 for the lagged variables, and we do not test
% for them here.
if any(isnan(data{enames{:}}(range).data(:)))
    error('Some variable values are missing in the database.')
end

% Set the number of exogenous variables.
xnbr = length(xnames);

% Test if we have a residual and get its name (-> rname).
if isequal(xnbr, 1)
    rname = M_.exo_names{strmatch(xnames{1}, M_.exo_names, 'exact')};
    irname = 1;
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

% Substitute parameters and variables.
for i=1:length(objNames)
    switch objTypes(i)
      case 1
        rhs = strrep(rhs, objNames{i}, sprintf('DynareModel.params(%u)', objIndex(i)));
      case {2,3}
        k = strmatch(objNames{i}, data.name, 'exact');
        if isempty(k)
            error('Variable %s is missing in the database.', objNames{i})
        end
        n = length(objNames{i});
        j = regexp(rhs, ['\<', objNames{i}, '\>']);
        if islaggedvariables
            jlag = regexp(rhs, ['\<(', objNames{i}, '\(-1\))\>']);
            if ~isempty(jlag)
                rhs = regexprep(rhs, ['\<(' objNames{i} '\(-1\))\>'], sprintf('data(1:end-1,%u)', k));
            end
            if ~isempty(setdiff(j, jlag))
                rhs = regexprep(rhs, ['\<' objNames{i} '\>'], sprintf('data(2:end,%u)', k));
            end
        else
            rhs = regexprep(rhs, ['\<' objNames{i} '\>'], sprintf('data(:,%u)', k));
        end
        if contains(lhs, objNames{i})
            if islaggedvariables
                lhs = strrep(lhs, objNames{i}, sprintf('data(2:end,%u)', k));
            else
                lhs = strrep(lhs, objNames{i}, sprintf('data(:,%u)', k));
            end
        end
    end
end

% Create a routine for evaluating the sum of squared residuals
ssrfun = ['ssr_' eqname];
fid = fopen([ssrfun '.m'], 'w');
fprintf(fid, 'function s = %s(params, data, DynareModel, DynareOutput)\n', ssrfun);
fprintf(fid, '\n');
fprintf(fid, '%% Evaluates the sum of square residuals for equation %s.\n', eqname);
fprintf(fid, '%% File created by Dynare (%s).\n', datestr(datetime));
fprintf(fid, '\n');
for i=1:length(ipnames_)
    fprintf(fid, 'DynareModel.params(%u) = params(%u);\n', ipnames_(i), i);
end
fprintf(fid, '\n');
fprintf(fid, 'DynareModel = pac.update.parameters(''%s'', DynareModel, DynareOutput);\n', ...
        pacmodl);
fprintf(fid, '\n');
fprintf(fid, 'r = %s-(%s);\n', lhs, rhs);
fprintf(fid, 's = .0;\n');
fprintf(fid, 'for i=1:%u\n', range.length()-islaggedvariables);
fprintf(fid, '    s = s + r(i)*r(i);\n');
fprintf(fid, 'end\n');
fclose(fid);

% Create a function handle returning the sum of square residuals for a given
% vector of parameters.
DATA = data(range).data;
ssr = @(p) feval(['ssr_' eqname], p, DATA, M_, oo_);

% Set initial condition.
params0 = cell2mat(struct2cell(params));

% Set optimization options.
options = optimset('Display','iter');

% Estimate the parameters by minimizing the sum of squared residuals.
pparams1 = fminunc(ssr, params0, options);

% Update M_.params
for i=1:length(pparams1)
    M_.params(ipnames_(i)) = pparams1(i);
end
