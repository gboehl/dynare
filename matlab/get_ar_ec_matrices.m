function get_ar_ec_matrices(var_model_name)
%function get_ar_ec_matrices(var_model_name)
%
% Returns the autoregressive and error correction matrices associated with the
% VAR specified by var_model_name. Output is stored in cellarray
% oo_.var.(var_model_name).ar, with oo_.var.(var_model_name).ar{i} being the
% AR matrix at time t-i (same holds for error correction matrices with ec
% replacing ar). Each AR (EC) matrix is stored with rows organized by the
% ordering of the equation tags found in M_.var.(var_model_name).eqtags and
% columns organized consistently.
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

%% Organize AR & EC matrices
assert(length(M_.var.(var_model_name).lhs) == rows(g1));

% Find RHS vars for AR & EC matrices
arRhsVars = [];
ecRhsVars = [];
lhs       = M_.var.(var_model_name).lhs;
for i = 1:length(M_.var.(var_model_name).rhs.vars_at_eq)
    vars = M_.var.(var_model_name).rhs.vars_at_eq{i}.var;
    rhsvars{i}.vars = [];
    rhsvars{i}.lags = [];
    rhsvars{i}.arRhsIdxs = [];
    rhsvars{i}.ecRhsIdxs = [];
    for j = 1:length(vars)
        if vars(j) <= M_.orig_endo_nbr
            % vars(j) is not an aux var
            if ismember(vars(j), lhs)
                arRhsVars = union(arRhsVars, vars(j), 'stable');
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs find(arRhsVars == vars(j))];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs -1];
            else
                ecRhsVars = union(ecRhsVars, vars(j), 'stable');
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs -1];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs find(ecRhsVars == vars(j))];
            end
        else
            % Search aux vars for matching lhs var
            lhsvaridx = findLhsInAuxVar(vars(j), lhs);
            if lhsvaridx >= 1
                arRhsVars = union(arRhsVars, lhsvaridx, 'stable');
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs find(arRhsVars == lhsvaridx)];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs -1];
            else
                % otherwise find endog that corresponds to this aux var
                varidx = findVarNoLag(vars(j));
                ecRhsVars = union(ecRhsVars, varidx, 'stable');
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs -1];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs find(ecRhsVars == varidx)];
            end
        end
    end
    rhsvars{i}.vars = vars;
    rhsvars{i}.lags = M_.var.(var_model_name).rhs.vars_at_eq{i}.lag;
end

% Initialize matrices
oo_.var.(var_model_name).ar = zeros(length(lhs), length(arRhsVars), M_.var.(var_model_name).max_lag);
oo_.var.(var_model_name).ec = zeros(length(lhs), length(ecRhsVars), M_.var.(var_model_name).max_lag);
oo_.var.(var_model_name).ar_idx = arRhsVars;
oo_.var.(var_model_name).ec_idx = ecRhsVars;
% Fill matrices
for i = 1:length(rhsvars)
    for j = 1:length(rhsvars{i}.vars)
        var = rhsvars{i}.vars(j);
        if rhsvars{i}.lags(j) == -1
            g1col = M_.lead_lag_incidence(1, var);
        else
            g1col = M_.lead_lag_incidence(2, var);
        end
        if g1col ~= 0 && any(g1(:, g1col))
            if rhsvars{i}.arRhsIdxs(j) > 0
                % Fill AR
                [lag, ndiffs] = findLagForVar(var, -rhsvars{i}.lags(j), 0, arRhsVars);
                if ndiffs >= 1
                    ndiffs = ndiffs - 1;
                end
                for k = 0:ndiffs
                    oo_.var.(var_model_name).ar(i, rhsvars{i}.arRhsIdxs(j), lag + k) = ...
                        oo_.var.(var_model_name).ar(i, rhsvars{i}.arRhsIdxs(j), lag + k) + (-1)^k * nchoosek(ndiffs,k) * g1(i, g1col);
                end
            elseif rhsvars{i}.ecRhsIdxs(j) > 0
                % Fill EC
                [lag, ndiffs] = findLagForVar(var, -rhsvars{i}.lags(j), 0, ecRhsVars);
                for k = 0:ndiffs
                    oo_.var.(var_model_name).ec(i, rhsvars{i}.ecRhsIdxs(j), lag + k) = ...
                        oo_.var.(var_model_name).ec(i, rhsvars{i}.ecRhsIdxs(j), lag + k) + (-1)^k * nchoosek(ndiffs,k) * g1(i, g1col);
                end
            else
                error('Shouldn''t arrive here');
            end
        end
    end
end
end



function lhsvaridx = findLhsInAuxVar(auxVar, lhsvars)

global M_

if auxVar <= M_.orig_endo_nbr
    lhsvaridx = -1;
    return
end

av = M_.aux_vars([M_.aux_vars.endo_index] == auxVar);

if ismember(av.orig_index, lhsvars)
    lhsvaridx = av.orig_index;
else
    lhsvaridx = findLhsInAuxVar(av.orig_index, lhsvars);
end
end



function idx = findVarNoLag(auxVar)

global M_

if auxVar <= M_.orig_endo_nbr
    error('Shouldn''t arrive here')
end

av = M_.aux_vars([M_.aux_vars.endo_index] == auxVar);

if ~isempty(av.unary_op_handle)
    idx = av.endo_index;
else
    if av.orig_index <= M_.orig_endo_nbr
        idx = av.orig_index;
    else
        idx = findVarNoLag(av.orig_index);
    end
end
end



function [lag, ndiffs] = findLagForVar(auxVar, lag, ndiffs, rhsVars)

global M_

if auxVar <= M_.orig_endo_nbr
    return
end

av = M_.aux_vars([M_.aux_vars.endo_index] == auxVar);

if av.type == 8
    ndiffs = ndiffs + 1;
end

if ismember(av.endo_index, rhsVars)
    if ~isempty(av.unary_op_handle) && (av.type == 8 || av.type == 9)
        lag = lag + abs(av.orig_lead_lag);
    end
elseif ismember(av.orig_index, rhsVars)
    if av.orig_index <= M_.orig_endo_nbr
        lag = lag + abs(av.orig_lead_lag);
    else
        [lag, ndiffs] = findLagForVar(av.orig_index, lag + 1, ndiffs, rhsVars);
    end
else
    if av.type == 8
        [lag, ndiffs] = findLagForVar(av.orig_index, lag, ndiffs, rhsVars);
    else
        [lag, ndiffs] = findLagForVar(av.orig_index, lag + 1, ndiffs, rhsVars);
    end
end
assert(lag > 0)
end

