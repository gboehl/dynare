function b = islagof(v1, v2)

% Returns true if and only if variable v1 is a lag of variable v2.
%
% INPUTS
% - v1        [string]     Variable name.
% - v2        [string]     Variable name.
%
% OUTPUTS
% - b         [logical]

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

global M_

if ~ischar(v1) && isscalar(v1) && isnumeric(v1) && isint(v1)
    id1 = v1;
    if id1>M_.endo_nbr
        error('First input must be an integer between 1 and %u', M_.endo_nbr)
    end
end

if ~ischar(v2) && isscalar(v2) && isnumeric(v2) && isint(v2)
    id2 = v2;
    if id2>M_.endo_nbr
        error('First input must be an integer between 1 and %u', M_.endo_nbr)
    end
end

% Set output default value
b = false;

% Get index of v1
if ischar(v1)
    id1 = find(strcmp(v1, M_.endo_names));
end

if isempty(id1)
    error('First input must be a variable name.')
end

% Get index of v2
if ischar(v2)
    id2 = find(strcmp(v2, M_.endo_names));
end

if isempty(id2)
    error('Second input must be a variable name.')
end

% A variable cannot be a lag of itself.
if isequal(id1, id2)
    return
end

% Are v1 and v2 auxiliary variables?
v1_is_an_auxiliary_variable = id1>M_.orig_endo_nbr;
v2_is_an_auxiliary_variable = id2>M_.orig_endo_nbr;

% At least one of the variables must be an auxiliary variable
if ~v1_is_an_auxiliary_variable && ~v2_is_an_auxiliary_variable
    return
end

% If v1 and v2 are auxiliary variables
if v1_is_an_auxiliary_variable && v2_is_an_auxiliary_variable
    auxinfo1 = M_.aux_vars(get_aux_variable_id(id1));
    auxinfo2 = M_.aux_vars(get_aux_variable_id(id2));
    isleadlag1 = ismember(auxinfo1.type, [1,0]);
    isleadlag2 = ismember(auxinfo2.type, [1,0]);
    % If v1 and v2 are lead/lag auxiliary variables
    if isleadlag1 && isleadlag2
        % If v1 and v2 are lead/lag of a common endogenous variable.
        if auxinfo1.orig_index && auxinfo2.orig_index
            if auxinfo1.orig_lead_lag<auxinfo2.orig_lead_lag
                b = true;
            end
        end
        return
    end
    isdiff1 = ismember(auxinfo1.type, [8,9]);
    isdiff2 = ismember(auxinfo1.type, [8,9]);
    if isdiff1 && isdiff2
        if isequal(auxinfo1.type, 9) && isequal(auxinfo2.type, 8)
            while isequal(auxinfo1.type, 9)
                if isequal(auxinfo1.orig_index, id2)
                    b = true;
                end
                auxinfo1 = M_.aux_vars(get_aux_variable_id(auxinfo1.orig_index));
            end
        end
        return
    end
end

if v1_is_an_auxiliary_variable && ~v2_is_an_auxiliary_variable
    auxinfo1 = M_.aux_vars(get_aux_variable_id(id1));
    % If v1 is not an auxiliary of type 1, it cannot be a lag of v2
    if isequal(auxinfo1.type, 1)
        if isequal(auxinfo1.orig_index, id2)
            b = true;
        end
    end
end

if ~v1_is_an_auxiliary_variable && v2_is_an_auxiliary_variable
    auxinfo2 = M_.aux_vars(get_aux_variable_id(id2));
    % If v2 is not an auxiliary of type 0, v1 cannot be a lag of v2
    if isequal(auxinfo2.type, 0)
        if isequal(auxinfo2.orig_index, id1)
            b = true;
        end
    end
end