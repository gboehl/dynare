function get_ar_matrices(var_model_name)
%function ar = get_ar_matrices(var_model_name)
% Gets the autoregressive matrices associated with the var specified by
% var_model_name. Output stored in cellarray oo_.var.(var_model_name).ar,
% with oo_.var.(var_model_name).ar(1) being the AR matrix at time t, 
% oo_.var.(var_model_name).ar(2) the AR matrix at time t-1, etc. Each
% AR matrix is stored with rows organized by the ordering of the equation
% tags found in M_.var.(var_model_name).eqtags and columns organized by
% M_.endo_names order
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
assert(ischar(var_model_name), 'The sole argument must be a string');
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
        continue
    end
end
g1 = -1 * g1(g1rows, :);

% Check for leads
if rows(M_.lead_lag_incidence) == 3
    idxs = M_.lead_lag_incidence(3, M_.lead_lag_incidence(3, :) ~= 0);
    assert(~any(g1(:, idxs)), ...
        ['You cannot have leads in the equations specified by ' strjoin(M_.var.(var_model_name).eqtags, ',')]);
end

%% Organize AR matrices
% Get LHS info
% NB: equations must have one endogenous variable on LHS
jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
lhs = getEquationsByTags(jsonmodel, 'name', M_.var.(var_model_name).eqtags);
lhsidxs = zeros(ntags, 1);
for i = 1:ntags
    idxs = strcmp(M_.endo_names, lhs{i});
    if any(idxs)
        lhsidxs(i) = find(idxs);
        continue
    end
end
assert(length(lhsidxs) == rows(g1));

% Initialize AR matrices
for i = 1:M_.max_endo_lag_orig+1
    oo_.var.(var_model_name).ar{i} = zeros(length(lhsidxs), M_.orig_endo_nbr);
    oo_.var.(var_model_name).artime{i} = 't';
    if i > 1
        oo_.var.(var_model_name).artime{i} = [oo_.var.(var_model_name).artime{i} '-' num2str(i-1)];
    end
end

for i = 1:2
    if i == 1
        baselag = 2;
    else
        baselag = 1;
    end
    for j = 1:size(M_.lead_lag_incidence, 2)
        if M_.lead_lag_incidence(i, j) ~= 0 && any(g1(:, M_.lead_lag_incidence(i, j)))
            if j > M_.orig_endo_nbr
                av = findauxvar(j);
                assert(~isempty(av));
                oo_.var.(var_model_name).ar{(av.orig_lead_lag * - 1) + baselag}(:, av.orig_index) = ...
                    g1(:, M_.lead_lag_incidence(i, j));
            else
                oo_.var.(var_model_name).ar{baselag}(:, j) = ...
                    g1(:, M_.lead_lag_incidence(i, j));
            end
        end
    end
end

sublhs = zeros(length(lhsidxs), M_.orig_endo_nbr);
for i = 1:length(lhsidxs)
    sublhs(i, lhsidxs(i)) = 1;
end
oo_.var.(var_model_name).ar{1} = oo_.var.(var_model_name).ar{1} + sublhs;
end

function out = findauxvar(endo_idx)
global M_
out = {};
for i=1:length(M_.aux_vars)
    if M_.aux_vars(i).endo_index == endo_idx
        out = M_.aux_vars(i);
        return
    end
end
end
