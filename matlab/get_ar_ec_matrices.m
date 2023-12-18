function get_ar_ec_matrices(model_name, model_type)
%function get_ar_ec_matrices(model_name, model_type)
%
% Returns the autoregressive matrix associated with the auxiliary model specified by
% model_name. Output is stored in cellarray oo_.(model_type).(model_name).ar,
% with oo_.(model_type).(model_name).ar(:,:,i) being the AR matrix at time t-i. Each
% AR matrix is stored with rows and columns organized by the ordering of the
% equation tags found in M_.(model_type).(model_name).eqtags.
% oo_.(model_type).(model_name).ec contains those entries that are not
% autoregressive.
%
% INPUTS
%
%   model_name   [string]        the name of the auxiliary model
%   model_type   [string]        the type of the auxiliary model ('var' or
%                                'trend_component'. If not passed, the
%                                value is set by the function; if a 'var'
%                                subfield is found, that is used. Otherwise
%                                'trend_component' is used (if it exists as
%                                a subfield of M_.
%
% OUTPUTS
%
%   NONE

% Copyright © 2018 Dynare Team
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

%% Check inputs
assert(nargin <= 2, 'This function requires one or two arguments');
assert(~isempty(model_name) && ischar(model_name), ...
    'The first argument must be a non-empty string');

if nargin < 2
    model_type = 'var';
    if ~(isfield(M_, model_type) && isfield(M_.(model_type), model_name))
        model_type = 'trend_component';
        if ~(isfield(M_, model_type) && isfield(M_.(model_type), model_name))
            error(['Could not find ' model_name ' in M_.var or ' ...
                'M_.trend_component. First declare it via the var_model ' ...
                'or trend_component_model statement.']);
        end
    end
else
    assert(~isempty(model_type) && ischar(model_type), ...
        'If provided, the second argument must be a non-empty string');
    if ~(isfield(M_, model_type) && isfield(M_.(model_type), model_name))
        error(['Could not find M_.' model_type '.' model_name ...
            '. First declare it via the var_model or ' ...
            'trend_component_model statement.']);
    end
end

%% Call Dynamic Function
[~, g1] = feval([M_.fname '.dynamic'], ...
    ones(max(max(M_.lead_lag_incidence)), 1), ...
    ones(1, M_.exo_nbr), ...
    M_.params, ...
    zeros(M_.endo_nbr, 1), ...
    1);

% Choose rows of Jacobian based on equation tags
ntags = length(M_.(model_type).(model_name).eqtags);
g1rows = zeros(ntags, 1);
for i = 1:ntags
    idxs = strcmp(M_.equations_tags(:, 3), M_.(model_type).(model_name).eqtags{i});
    if any(idxs)
        g1rows(i) = M_.equations_tags{idxs, 1};
    end
end
g1 = -1 * g1(g1rows, :);

% Check for leads
if rows(M_.lead_lag_incidence) == 3
    idxs = M_.lead_lag_incidence(3, M_.lead_lag_incidence(3, :) ~= 0);
    assert(~any(any(g1(g1rows, idxs))), ...
        ['You cannot have leads in the equations specified by ' strjoin(M_.(model_type).(model_name).eqtags, ',')]);
end

%% Organize AR & EC matrices
assert(length(M_.(model_type).(model_name).lhs) == rows(g1));

% Find RHS vars for AR & EC matrices
ecRhsVars = [];
lhs       = M_.(model_type).(model_name).lhs;
rhsvars   = cell(length(lhs), 1);
for i = 1:length(M_.(model_type).(model_name).rhs.vars_at_eq)
    vars = M_.(model_type).(model_name).rhs.vars_at_eq{i}.var;
    rhsvars{i}.vars = vars;
    rhsvars{i}.lags = M_.(model_type).(model_name).rhs.vars_at_eq{i}.lag;
    rhsvars{i}.arRhsIdxs = [];
    rhsvars{i}.ecRhsIdxs = [];
    rhsvars{i}.ecRhsVars = [];
    for j = 1:length(vars)
        if vars(j) <= M_.orig_endo_nbr
            % vars(j) is not an aux var
            if ismember(vars(j), lhs)
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs find(lhs == vars(j))];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs -1];
                rhsvars{i}.ecRhsVars = [rhsvars{i}.ecRhsVars -1];
            else
                ecRhsVars = union(ecRhsVars, vars(j));
                rhsvars{i}.ecRhsVars = [rhsvars{i}.ecRhsVars vars(j)];
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs -1];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs find(rhsvars{i}.ecRhsVars == vars(j))];
            end
        else
            % Search aux vars for matching lhs var
            lhsvaridx = findLhsInAuxVar(vars(j), lhs);
            if lhsvaridx >= 1
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs find(lhs == lhsvaridx)];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs -1];
                rhsvars{i}.ecRhsVars = [rhsvars{i}.ecRhsVars -1];
            else
                % otherwise find endog that corresponds to this aux var
                varidx = findVarNoLag(vars(j));
                ecRhsVars = union(ecRhsVars, varidx);
                rhsvars{i}.ecRhsVars = [rhsvars{i}.ecRhsVars varidx];
                rhsvars{i}.arRhsIdxs = [rhsvars{i}.arRhsIdxs -1];
                rhsvars{i}.ecRhsIdxs = [rhsvars{i}.ecRhsIdxs find(rhsvars{i}.ecRhsVars == varidx)];
            end
        end
    end
end

[rhsvars, ecRhsVars] = reorderECvars(rhsvars, ecRhsVars, lhs);

% Initialize matrices
oo_.(model_type).(model_name).ar = zeros(length(lhs), length(lhs), max(M_.(model_type).(model_name).max_lag));
oo_.(model_type).(model_name).ec = zeros(length(lhs), length(ecRhsVars), 1);
oo_.(model_type).(model_name).ar_idx = lhs;
oo_.(model_type).(model_name).ec_idx = ecRhsVars;

% Fill matrices
for i = 1:length(lhs)
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
                lag = findLagForVar(var, -rhsvars{i}.lags(j), 0, lhs);
                oo_.(model_type).(model_name).ar(i, rhsvars{i}.arRhsIdxs(j), lag) = ...
                    oo_.(model_type).(model_name).ar(i, rhsvars{i}.arRhsIdxs(j), lag) + g1(i, g1col);
            elseif rhsvars{i}.ecRhsIdxs(j) > 0
                % Fill EC
                lag = findLagForVar(var, -rhsvars{i}.lags(j), 0, ecRhsVars);
                if lag==1
                    if size(oo_.(model_type).(model_name).ec, 3) < lag
                        oo_.(model_type).(model_name).ec(i, rhsvars{i}.ecRhsIdxs(j), lag) = 0;
                    end
                    oo_.(model_type).(model_name).ec(i, rhsvars{i}.ecRhsIdxs(j), lag) = ...
                        oo_.(model_type).(model_name).ec(i, rhsvars{i}.ecRhsIdxs(j), lag) + g1(i, g1col);
                end
            else
                error('Shouldn''t arrive here');
            end
        end
    end
end
end


function [rhsvars, ecRhsVarsReordered] = reorderECvars(rhsvars, ecRhsVars, lhs)

global M_

ecRhsVarsReordered = [];
for i = 1:length(lhs)
    av = M_.aux_vars([M_.aux_vars.endo_index] == lhs(i));
    if ~isempty(av)
        var = ecRhsVars(ecRhsVars == av.orig_index);
        if ~isempty(var)
            ecRhsVarsReordered = [ecRhsVarsReordered var];
             for j = 1:length(rhsvars)
                 rhsidx = find(rhsvars{j}.ecRhsVars == var);
                 if ~isempty(rhsidx)
                     rhsvars{j}.ecRhsIdxs(rhsidx) = length(ecRhsVarsReordered);
                 end
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

if av.type == 8 || av.type == 10
    if ismember(av.endo_index, lhsvars)
        lhsvaridx = av.endo_index;
    else
        lhsvaridx = findLhsInAuxVar(av.orig_index, lhsvars);
    end
else
    if ismember(av.orig_index, lhsvars)
        lhsvaridx = av.orig_index;
    else
        lhsvaridx = findLhsInAuxVar(av.orig_index, lhsvars);
    end
end
end



function idx = findVarNoLag(auxVar)

global M_

if auxVar <= M_.orig_endo_nbr
    error('Shouldn''t arrive here')
end

av = M_.aux_vars([M_.aux_vars.endo_index] == auxVar);

if av.type == 8 || av.type == 10
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
    assert(lag > 0)
    return
end

av = M_.aux_vars([M_.aux_vars.endo_index] == auxVar);

if av.type == 8
    ndiffs = ndiffs + 1;
end

if ismember(av.endo_index, rhsVars)
    if av.type == 8 || av.type == 9
        lag = lag + abs(av.orig_lead_lag);
    end
elseif ismember(av.orig_index, rhsVars)
    if av.orig_index <= M_.orig_endo_nbr
        lag = lag + abs(av.orig_lead_lag);
    else
        [lag, ndiffs] = findLagForVar(av.orig_index, lag + max(1, abs(av.orig_lead_lag)), ndiffs, rhsVars);
    end
else
    if av.type == 8
        [lag, ndiffs] = findLagForVar(av.orig_index, lag, ndiffs, rhsVars);
    else
        [lag, ndiffs] = findLagForVar(av.orig_index, lag + max(1, abs(av.orig_lead_lag)), ndiffs, rhsVars);
    end
end
assert(lag > 0)
end

