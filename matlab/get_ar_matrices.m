function get_ar_matrices(var_model_name)

% Gets the autoregressive matrices associated with the var specified by var_model_name.
% Output stored in cellarray oo_.var.(var_model_name).AutoregressiveMatrices, with
% oo_.var.(var_model_name).AutoregressiveMatrices{i} being the AR matrix at time t-i. Each
% AR matrix is stored with rows organized by the ordering of the equation tags
% found in M_.var.(var_model_name).eqtags and columns organized consistently.
%
% INPUTS
%
%   var_model_name   [string]        the name of the VAR model
%
% OUTPUTS
%
%   NONE

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

%% Check inputs and initialize output
assert(nargin == 1, 'This function requires one argument');
assert(~isempty(var_model_name) && ischar(var_model_name), ...
    'The sole argument must be a non-empty string');
if ~isfield(M_.var, var_model_name)
    error(['Could not find ' var_model_name ' in M_.var. ' ...
        'First declare it via the var_model statement.']);
end

%% Call Dynamic Function
[junk, g1] = feval([M_.fname '_dynamic'], ...
    ones(max(max(M_.lead_lag_incidence)), 1), ...
    ones(1, M_.exo_nbr), ...
    M_.params, ...
    zeros(M_.endo_nbr, 1), ...
    1);

% Choose rows of Jacobian based on equation tags
ntags = length(M_.var.(var_model_name).eqtags);
g1rows = zeros(ntags, 1);
for i = 1:ntags
    idxs = strcmp(M_.equations_tags(:, 3), M_.var.(var_model_name).eqtags{i});
    if any(idxs)
        g1rows(i) = M_.equations_tags{idxs, 1};
    end
end
g1 = -1 * g1(g1rows, :);

% Check for leads
if rows(M_.lead_lag_incidence) == 3
    idxs = M_.lead_lag_incidence(3, M_.lead_lag_incidence(3, :) ~= 0);
    assert(~any(any(g1(g1rows, idxs))), ...
        ['You cannot have leads in the equations specified by ' strjoin(M_.var.(var_model_name).eqtags, ',')]);
end

%% Organize AR matrices
assert(length(M_.var.(var_model_name).lhs) == rows(g1));

% Initialize AR matrices
% ECM
rhsvars = [];
rhslag = [];
maxlag = max(M_.var.(var_model_name).rhs.lag);
for i = 1:length(M_.var.(var_model_name).rhs.vars_at_eq)
    rhsvars = union(rhsvars, M_.var.(var_model_name).rhs.vars_at_eq{i}.var);
    rhslag = union(rhslag, M_.var.(var_model_name).rhs.vars_at_eq{i}.lag);
end
if ~isempty(rhslag)
    maxlag = max(maxlag, max(abs(rhslag)));
end

Bvars = setdiff(rhsvars, M_.var.(var_model_name).lhs);
orig_diff_var_vec = M_.var.(var_model_name).orig_diff_var(M_.var.(var_model_name).diff);
diff_vars = M_.var.(var_model_name).lhs(M_.var.(var_model_name).lhs > M_.orig_endo_nbr);
drop_avs_related_to = [];
for i = 1:length(diff_vars)
    av = M_.aux_vars([M_.aux_vars.endo_index] == diff_vars(i));
    assert(any(orig_diff_var_vec == av.orig_index));
    if av.type == 8
        drop_diff_avs_related_to = [drop_avs_related_to av.orig_index];
    end
end
keep = true(length(Bvars), 1);
Bvars_diff_index = zeros(length(Bvars), 1);
for i = sum(Bvars <= M_.orig_endo_nbr)+1:length(Bvars)
    av = M_.aux_vars([M_.aux_vars.endo_index] == Bvars(i));
    if av.type == 8
        if any(drop_diff_avs_related_to == av.orig_index)
            assert(any(orig_diff_var_vec == av.orig_index));
            keep(i) = false;
        else
            if any(Bvars_diff_index == av.orig_index)
                keep(i) = false;
            else
                Bvars_diff_index(i) = av.orig_index;
            end
        end
    end
end
Bvars = Bvars(keep);
Bvars_diff_index = Bvars_diff_index(keep);
oo_.var.(var_model_name).ecm_idx = Bvars;

% AR
narvars = length(M_.var.(var_model_name).lhs);
nothvars = length(Bvars);
for i = 1:maxlag + 1
    oo_.var.(var_model_name).AutoregressiveMatrices{i} = zeros(narvars, narvars);
    oo_.var.(var_model_name).ecm{i} = zeros(narvars, nothvars);
end

ecm_assigned = false;
for i = 1:2
    if i == 1
        baselag = 2;
    else
        baselag = 1;
    end
    for j = 1:size(M_.lead_lag_incidence, 2)
        if M_.lead_lag_incidence(i, j) ~= 0 && any(g1(:, M_.lead_lag_incidence(i, j)))
            if j > M_.orig_endo_nbr
                av = M_.aux_vars([M_.aux_vars.endo_index] == j);
                assert(~isempty(av));
                if av.type == 8
                    col = M_.var.(var_model_name).orig_diff_var == av.orig_index;
                    if any(col)
                        assert(any(orig_diff_var_vec == av.orig_index));
                        oo_.var.(var_model_name).AutoregressiveMatrices{(av.orig_lead_lag * - 1) + baselag}(:, col) = ...
                            g1(:, M_.lead_lag_incidence(i, j));
                    else
                        col = Bvars_diff_index == av.orig_index;
                        ecm_assigned = true;
                        oo_.var.(var_model_name).ecm{(av.orig_lead_lag * - 1) + baselag}(:, col) = ...
                            g1(:, M_.lead_lag_incidence(i, j));
                    end
                else
                    col = M_.var.(var_model_name).lhs == av.orig_index;
                    if any(col)
                        oo_.var.(var_model_name).AutoregressiveMatrices{(av.orig_lead_lag * - 1) + baselag}(:, col) = ...
                            g1(:, M_.lead_lag_incidence(i, j));
                    else
                        col = Bvars == av.orig_index;
                        ecm_assigned = true;
                        oo_.var.(var_model_name).ecm{(av.orig_lead_lag * - 1) + baselag}(:, col) = ...
                            g1(:, M_.lead_lag_incidence(i, j));
                    end
                end
            else
                col = M_.var.(var_model_name).lhs == j;
                if ~any(col)
                    col = Bvars == j;
                    ecm_assigned = true;
                    oo_.var.(var_model_name).ecm{baselag}(:, col) = ...
                        g1(:, M_.lead_lag_incidence(i, j));
                else
                    oo_.var.(var_model_name).AutoregressiveMatrices{baselag}(:, col) = ...
                        g1(:, M_.lead_lag_incidence(i, j));
                end
            end
        end
    end
end

for i = 1:length(M_.var.(var_model_name).lhs)
    oo_.var.(var_model_name).AutoregressiveMatrices{1}(i, M_.var.(var_model_name).lhs == M_.var.(var_model_name).lhs(i)) = ...
        oo_.var.(var_model_name).AutoregressiveMatrices{1}(i, M_.var.(var_model_name).lhs == M_.var.(var_model_name).lhs(i)) + 1;
end

if any(oo_.var.(var_model_name).AutoregressiveMatrices{1}(:))
    error('This is not a VAR model! Contemporaneous endogenous variables are not allowed.')
end

% Remove time t matrix for autoregressive part
oo_.var.(var_model_name).AutoregressiveMatrices = oo_.var.(var_model_name).AutoregressiveMatrices(2:end);

% Remove error correction matrices if never assigned
if ~ecm_assigned
    oo_.var.(var_model_name) = rmfield(oo_.var.(var_model_name), 'ecm');
    oo_.var.(var_model_name) = rmfield(oo_.var.(var_model_name), 'ecm_idx');
end
end
