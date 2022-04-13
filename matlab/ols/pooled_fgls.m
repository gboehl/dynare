function pooled_fgls(ds, param_common, param_regex, eqtags, model_name, param_names, ds_range)
% function pooled_fgls(ds, param_common, param_regex, eqtags, model_name, param_names, ds_range)
% Run Pooled FGLS
%
% INPUTS
%   ds            [dseries]      data to use in estimation
%   param_common  [cellstr]      List of values to insert into param_regex,
%                                e.g. country codes {'FR', 'DE', 'IT'}
%   param_regex   [cellstr]      Where '*' should be replaced by the first
%                                value in param_common
%   eqtags        [cellstr]      names of equation tags to estimate. If empty,
%                                estimate all equations
%   model_name    [string]       name to use in oo_ and inc file
%
%   param_names   [cellstr]      list of parameters to estimate (if
%                                empty, estimate all) (may contain regex
%                                to match param_regex)
%   ds_range      [dates]   range of dates to use in estimation
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright © 2017-2019 Dynare Team
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

global M_ oo_

%% Check input arguments
if nargin < 1 || nargin > 7
    error('Incorrect number of arguments')
end

if isempty(ds) || ~isdseries(ds)
    error('The first argument must be a dseries');
end

if nargin < 7
    ds_range = ds.dates;
end

if nargin < 6
    param_names = {};
end

if nargin < 5
    model_name = '';
end

if nargin < 4
    eqtags = {};
end

maxit = 100;
tol = 1e-6;

%% Common work between pooled_ols and pooled_fgls
[Y, X, pbeta, residnames, country_name, model_name] = pooled_ols(ds, param_common, param_regex, true, eqtags, model_name, param_names, ds_range);

%% Estimation
neqs = length(residnames);
oo_.pooled_fgls.(model_name).dof = size(X,1)/neqs;
beta0 = oo_.pooled_fgls.(model_name).beta;
for i = 1:maxit
    resid = Y - X * beta0;
    resid = reshape(resid, oo_.pooled_fgls.(model_name).dof, neqs);
    vcv = resid'*resid/oo_.pooled_fgls.(model_name).dof;
    kLeye = kron(inv(chol(vcv))', eye(oo_.pooled_fgls.(model_name).dof));
    [q, r] = qr(kLeye*X, 0);
    oo_.pooled_fgls.(model_name).beta = r\(q'*kLeye*Y);
    if max(abs(beta0 - oo_.pooled_fgls.(model_name).beta)) < tol
        break
    end
    beta0 = oo_.pooled_fgls.(model_name).beta;
    if i == maxit
        warning('maximum nuber of iterations reached')
    end
end

% Set appropriate entries in M_.Sigma_e
idxs = zeros(neqs, 1);
for i = 1:neqs
    idxs(i) = find(strcmp(residnames{i}, M_.exo_names));
end
M_.Sigma_e(idxs, idxs) = vcv;

regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
assigned_idxs = false(size(pbeta));
incidxs = [];
for i = 1:length(param_regex)
    beta_idx = strcmp(pbeta, strrep(param_regex{i}, '*', country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.pooled_fgls.(model_name).beta(beta_idx);
    if isempty(eqtags) && isempty(param_names)
        assert(~isempty(value));
    end
    if ~isempty(value)
        idxs = find(~cellfun(@isempty, regexp(M_.param_names, strrep(param_regex{i}, '*', regexcountries))))';
        incidxs = [incidxs idxs];
        M_.params(idxs) = value;
    end
end
idxs = find(assigned_idxs == 0);
values = oo_.pooled_fgls.(model_name).beta(idxs);
names = pbeta(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    incidxs = [incidxs find(strcmp(M_.param_names, names{i}))];
    M_.params(incidxs(end)) = values(i);
end

% Write .inc file
write_param_init_inc_file('pooled_fgls', model_name, incidxs, M_.params(incidxs));

end
