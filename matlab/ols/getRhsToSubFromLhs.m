function lhssub = getRhsToSubFromLhs(ds, rhs, regex, splits, pnames)
%function lhssub = getRhsToSubFromLhs(ds, rhs, regex, splits, pnames)
% Helper function that identifies variables on RHS that need to be
% subtracted from LHS of OLS-style equation
%
% INPUTS
%   ds                [dseries]     data
%   rhs               [string]      RHS as a string
%   regex             [string]      regex expressing valid list of variables
%   splits            [cell string] strings to split out of equation on RHS
%   pnames            [cell string] parameter names
%
% OUTPUTS
%   lhssub            [dseries]    summed data to subtract from LHS
%
% SPECIAL REQUIREMENTS
%   none

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

global M_

assert(isdseries(ds), 'The first argument must be a dseries');
assert(ischar(rhs), 'The second argument must be a string');
assert(ischar(regex), 'The third argument must be a string');
assert(iscellstr(splits), 'The fourth argument must be a cell');
assert(iscellstr(splits), 'The fourth argument must be a cell');

lhssub = dseries();
rhs_ = strsplit(rhs, splits);
for j = 1:length(rhs_)
    rhsj = rhs_{j};
    while ~isempty(rhsj)
        minusstr = '';
        if strcmp(rhsj(1), '-') || strcmp(rhsj(1), '+')
            if length(rhsj) == 1
                break
            end
            if strcmp(rhsj(1), '-')
                minusstr = '-';
            end
            rhsj = rhsj(2:end);
        end
        str = getStrMoveRight(rhsj);
        if ~isempty(str)
            try
                lhssub = lhssub + eval(regexprep([minusstr str], regex, 'ds.$&'));
            catch
                if ~any(strcmp(M_.exo_names, str)) && ~any(strcmp(pnames, str))
                    error(['getRhsToSubFromLhs: problem evaluating ' minusstr str]);
                end
            end
            rhsj = rhsj(length(str)+1:end);
        end
    end
end
if ~isempty(lhssub)
    assert(lhssub.vobs == 1, 'error in getRhsToSubFromLhs');
    lhssub.rename_(lhssub.name{:}, 'summed_rhs');
end
end