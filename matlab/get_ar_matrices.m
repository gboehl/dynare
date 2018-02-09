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
    assert(~any(g1(:, idxs)), ...
        ['You cannot have leads in the equations specified by ' strjoin(M_.var.(var_model_name).eqtags, ',')]);
end

%% Organize AR matrices
assert(length(M_.var.(var_model_name).lhs) == rows(g1));

% Initialize AR matrices
for i = 1:max(M_.var.(var_model_name).rhs.lag)+1
    oo_.var.(var_model_name).ar{i} = zeros(length(M_.var.(var_model_name).lhs), length(M_.var.(var_model_name).lhs));
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
                av = M_.aux_vars([M_.aux_vars.endo_index] == j);
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

for i = 1:length(M_.var.(var_model_name).lhs)
    oo_.var.(var_model_name).ar{1}(i, M_.var.(var_model_name).lhs(i)) = oo_.var.(var_model_name).ar{1}(i, M_.var.(var_model_name).lhs(i)) + 1;
end

if any(oo_.var.(var_model_name).ar{1}(:))
    error('This is not a VAR model! Contemporaneous endogenous variables are not allowed.')
end

% Remove first matrix
oo_.var.(var_model_name).ar = oo_.var.(var_model_name).ar(2:end);