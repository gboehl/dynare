function ds = dyn_ols(ds, fitted_names_dict, eqtags, model_name, param_names, ds_range, update_params)
% function varargout = dyn_ols(ds, fitted_names_dict, eqtags, model_name, param_names, ds_range)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds                [dseries]         data
%   fitted_names_dict [cell]            Nx2 or Nx3 cell array to be used in naming fitted
%                                       values; first column is the equation tag,
%                                       second column is the name of the
%                                       associated fitted value, third column
%                                       (if it exists) is the function name of
%                                       the transformation to perform on the
%                                       fitted value.
%   eqtags            [cellstr]         names of equation tags to estimate. If empty,
%                                       estimate all equations
%   model_name        [celltsr]         name to use in oo_ and inc file (must be
%                                       same size as eqtags)
%   param_names       [cell of cellstr] list of parameters to estimate by eqtag
%                                       (if empty, estimate all)
%   ds_range          [dates]           range of dates to use in estimation
%   update_params     [logical]         If false M_.params will not be estimation.
%                                       If true (default) M_.params will be updated.
%
% OUTPUTS
%   ds                [dseries]    data updated with fitted values
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright (C) 2017-2021 Dynare Team
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

global M_ oo_ options_

if nargin < 1 || nargin > 7
    error('dyn_ols() takes between 1 and 7 arguments')
end

if isempty(ds) || ~isdseries(ds)
    error('dyn_ols: the first argument must be a dseries')
end

if nargin == 6 && ~isempty(ds_range)
    update_params = true;
elseif nargin < 6
    ds_range = ds.dates;
    update_params = true;
else
    if isempty(ds_range)
        ds_range = ds.dates;
    else
        if ds_range(1) < ds.firstdate || ds_range(end) > lastdate(ds)
            error('There is a problem with the 6th argument: the date range does not correspond to that of the dseries')
        end
    end
    update_params = true;
end

if nargin < 5
    param_names = {};
else
    if ~isempty(param_names)
        if ~iscell(param_names) || (~isempty(eqtags) && length(param_names) ~= length(eqtags))
            error('The 5th argument, if provided, must be a cell of the same length as the eqtags argument')
        end
        for i = 1:length(param_names)
            if ~iscellstr(param_names{i})
                error('every entry of param_names must be a cellstr')
            end
        end
    end
end

if nargin < 3
    eqtags = {};
end

if nargin < 2
    fitted_names_dict = {};
else
    assert(isempty(fitted_names_dict) || ...
        (iscell(fitted_names_dict) && ...
        (size(fitted_names_dict, 2) == 2 || size(fitted_names_dict, 2) == 3)), ...
        'dyn_ols: the second argument must be an Nx2 or Nx3 cell array');
end

%% Get Equation(s)
ast = get_ast(eqtags);

%% Set model_name
if nargin < 4
    model_name = cell(length(ast), 1);
else
    if isempty(model_name)
        model_name = repmat({''}, length(ast), 1);
    else
        if ~iscellstr(model_name) || length(model_name) ~= length(ast)
            error('The length of the 4th argument must be a cellstr with length equal to the number of equations estimated')
        end
        for i = 1:length(model_name)
            if ~isvarname(model_name{i})
                error('Every entry in the 4th argument must be a valid string');
            end
        end
    end
end

%% Parse equations
[Y, lhssub, X, fp, lp] = common_parsing(ds(ds_range), ast, true, param_names);

%% Loop over equations
for i = 1:length(Y)
    pnames = X{i}.name;
    [nobs, nvars] = size(X{i}.data);

    if ~isempty(model_name{i})
        tag = model_name{i};
    else
        if isfield(ast{i}, 'tags') && isfield(ast{i}.tags, 'name')
            tag = ast{i}.tags.('name');
        else
            tag = ['eq_line_no_' num2str(ast{i}.line)];
        end
    end

    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    oo_.ols.(tag).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X{i}.data, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tag).beta = r\(q'*Y{i}.data);
    oo_.ols.(tag).param_idxs = zeros(length(pnames), 1);
    for j = 1:length(pnames)
        if ~strcmp(pnames{j}, 'intercept')
            oo_.ols.(tag).param_idxs(j) = find(strcmp(M_.param_names, pnames{j}));
            if update_params
                M_.params(oo_.ols.(tag).param_idxs(j)) = oo_.ols.(tag).beta(j);
            end
        end
    end

    % Write .inc file
    write_param_init_inc_file('ols', tag, oo_.ols.(tag).param_idxs, oo_.ols.(tag).beta);

    % Yhat
    idx = 0;
    yhatname = [tag '_FIT'];
    if ~isempty(fitted_names_dict)
        idx = strcmp(fitted_names_dict(:,1), tag);
        if any(idx)
            yhatname = fitted_names_dict{idx, 2};
        end
    end
    oo_.ols.(tag).Yhat = dseries(X{i}.data*oo_.ols.(tag).beta, fp{i}, yhatname);

    % Residuals
    oo_.ols.(tag).resid = Y{i} - oo_.ols.(tag).Yhat;

    % Correct Yhat reported back to user
    Y{i} = Y{i} + lhssub{i};
    oo_.ols.(tag).Yobs = Y{i};
    oo_.ols.(tag).Yhat = oo_.ols.(tag).Yhat + lhssub{i};
    oo_.ols.(tag).YhatOrig = oo_.ols.(tag).Yhat;

    % Apply correcting function for Yhat if it was passed
    if any(idx) ...
            && length(fitted_names_dict(idx, :)) == 3 ...
            && ~isempty(fitted_names_dict{idx, 3})
        oo_.ols.(tag).Yhat = ...
            feval(fitted_names_dict{idx, 3}, oo_.ols.(tag).Yhat);
    end
    ds.(oo_.ols.(tag).Yhat.name{:}) = oo_.ols.(tag).Yhat;

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tag).resid.data'*oo_.ols.(tag).resid.data;
    oo_.ols.(tag).s2 = SS_res/oo_.ols.(tag).dof;

    % R^2
    ym = Y{i}.data - mean(Y{i});
    SS_tot = ym'*ym;
    oo_.ols.(tag).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tag).adjR2 = oo_.ols.(tag).R2 - (1 - oo_.ols.(tag).R2)*(nvars-1)/(oo_.ols.(tag).dof);

    % Durbin-Watson
    ediff = oo_.ols.(tag).resid.data(2:nobs) - oo_.ols.(tag).resid.data(1:nobs-1);
    oo_.ols.(tag).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tag).stderr = sqrt(oo_.ols.(tag).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tag).tstat = oo_.ols.(tag).beta./oo_.ols.(tag).stderr;

    %% Print Output
    if ~options_.noprint
        if nargin == 3
            title = ['OLS Estimation of equation ''' tag ''' [name = ''' tag ''']'];
        else
            title = ['OLS Estimation of equation ''' tag ''''];
        end

        preamble = {['Dependent Variable: ' Y{i}.name{:}], ...
            sprintf('No. Independent Variables: %d', nvars), ...
            sprintf('Observations: %d from %s to %s\n', nobs, fp{i}.char, lp{i}.char)};

        afterward = {sprintf('R^2: %f', oo_.ols.(tag).R2), ...
            sprintf('R^2 Adjusted: %f', oo_.ols.(tag).adjR2), ...
            sprintf('s^2: %f', oo_.ols.(tag).s2), ...
            sprintf('Durbin-Watson: %f', oo_.ols.(tag).dw)};

        dyn_table(title, preamble, afterward, pnames, ...
            {'Estimates','t-statistic','Std. Error'}, 4, ...
            [oo_.ols.(tag).beta oo_.ols.(tag).tstat oo_.ols.(tag).stderr]);
    end
end
end
