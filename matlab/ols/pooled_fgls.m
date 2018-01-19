function pooled_fgls(ds, param_common, param_regex, eqtags)
% function pooled_fgls(ds, param_common, param_regex, eqtags)
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
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must be run with the option: json=parse

% Copyright (C) 2017-2018 Dynare Team
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

if nargin < 4
    pooled_ols(ds, param_common, param_regex, true);
else
    param_regex = pooled_ols(ds, param_common, param_regex, true, eqtags);
end

oo_.sur.dof = length(oo_.pooled_fgls.sample_range);
resid = oo_.pooled_fgls.Y - oo_.pooled_fgls.X * oo_.pooled_fgls.beta;
resid = reshape(resid, oo_.sur.dof, length(oo_.pooled_fgls.residnames));
M_.Sigma_e = resid'*resid/oo_.sur.dof;

kLeye = kron(chol(inv(M_.Sigma_e)), eye(oo_.sur.dof));
[q, r] = qr(kLeye*oo_.pooled_fgls.X, 0);
oo_.pooled_fgls.beta = r\(q'*kLeye*oo_.pooled_fgls.Y);

regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
pbeta = oo_.pooled_fgls.pbeta;
assigned_idxs = false(size(pbeta));
for i = 1:length(param_regex)
    beta_idx = strcmp(pbeta, strrep(param_regex{i}, '*', oo_.pooled_fgls.country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.pooled_fgls.beta(beta_idx);
    if isempty(eqtags)
        assert(~isempty(value));
    end
    if ~isempty(value)
        M_.params(~cellfun(@isempty, regexp(M_.param_names, ...
            strrep(param_regex{i}, '*', regexcountries)))) = value;
    end
end
idxs = find(assigned_idxs == 0);
values = oo_.pooled_fgls.beta(idxs);
names = pbeta(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    M_.params(strcmp(M_.param_names, names{i})) = values(i);
end

oo_.pooled_fgls = rmfield(oo_.pooled_fgls, 'X');
oo_.pooled_fgls = rmfield(oo_.pooled_fgls, 'Y');
oo_.pooled_fgls = rmfield(oo_.pooled_fgls, 'residnames');
oo_.pooled_fgls = rmfield(oo_.pooled_fgls, 'pbeta');
oo_.pooled_fgls = rmfield(oo_.pooled_fgls, 'country_name');

end
