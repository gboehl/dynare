function retval = getStrMoveLeft(str)
%function retval = getStrMoveLeft(str)
% Helper function common to OLS routines. Given a string finds expr
% moving from right to left (potentially contained in parenthesis) until
% it gets to the + or -. e.g., given: str = 
% 'res_9+DE_SIN*param+log(DE_SIN(-1))', returns 'log(DE_SIN)'
%
% INPUTS
%   str       [string]    string
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017 Dynare Team
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

mathops = '[\+\*\^\-\/]';

if str(end) ~= ')'
    mathidxs = regexp(str, mathops);
    if isempty(mathidxs)
        retval = str;
    else
        retval = str(max(mathidxs)+1:end);
    end
else
    closedidxs = strfind(str, ')');
    closedidxs = [(length(closedidxs):-1:1)' closedidxs'];
    openidxs = strfind(str, '(');
    openidxs = [(length(openidxs):-1:1)' openidxs'];
    assert(rows(closedidxs) == rows(openidxs));
    for i = rows(openidxs):-1:1
        openparenidx = find(openidxs(i, 2) < closedidxs(:, 2), 1, 'first');
        if openidxs(i, 1) == closedidxs(openparenidx, 1)
            break;
        end
    end
    retval = str(max(regexp(str(1:openidxs(openparenidx,2)), mathops))+1:closedidxs(end));
end
