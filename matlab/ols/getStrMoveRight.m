function retval = getStrMoveRight(str)
%function retval = getStrMoveRight(str)
% Helper function common to OLS routines. Given a string finds expr
% moving from left to right (potentially contained in parenthesis) until
% it gets to the + or -. e.g., given: str = 
% '(log(DE_SIN(-1))-log(DE_SIN))+anotherparam1*log(DE_SIN)', returns
% '(log(DE_SIN(-1))-log(DE_SIN))'
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

if isempty(str)
    retval = '';
    return
end

mathidxs = regexp(str, '[\+\-]');
openidxs = strfind(str, '(');
if isempty(openidxs) ...
        || (~isempty(mathidxs) ...
        && min(mathidxs) < min(openidxs))
    if isempty(mathidxs)
        retval = str;
    else
        retval = str(1:min(mathidxs)-1);
    end
else
    openidxs = [(1:length(openidxs))' openidxs'];
    closedidxs = strfind(str, ')');
    closedidxs = [(1:length(closedidxs))' closedidxs'];
    assert(rows(openidxs) == rows(closedidxs));
    for i = 1:rows(closedidxs)
        closedparenidx = sum(openidxs(:, 2) < closedidxs(i, 2));
        if openidxs(closedparenidx, 1) == closedidxs(i, 1)
            break;
        end
    end
    retval = str(1:closedidxs(closedparenidx, 2));
    if length(str) > closedidxs(closedparenidx, 2) + 1
        if any(strfind('*^/', str(closedidxs(closedparenidx, 2) + 1)))
            retval = [retval ...
                str(closedidxs(closedparenidx, 2) + 1) ...
                getStrMoveRight(str(closedidxs(closedparenidx, 2) + 2:end))];
        end
    end
end
