function p = subsref(o, S)

% Overload subsref method.

% Copyright © 2023 Dynare Team
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

switch S(1).type
  case '.'
    if ismember(S(1).subs, {'p1','p2','p3','p4','p5','p6','p7','lb','ub'})
        p = builtin('subsref', o, S(1));
    elseif ismember(S(1).subs, {'draw','length'})
        p = feval(S(1).subs, o);
    elseif ismember(S(1).subs, {'draws', 'density', 'densities', 'moments', 'admissible'})
        p = feval(S(1).subs, o , S(2).subs{:});
    elseif ismember(S(1).subs, {'mean', 'median', 'variance', 'mode'})
        if (length(S)==2 && isempty(S(2).subs)) || length(S)==1
            p = feval(S(1).subs, o);
        else
            p = feval(S(1).subs, o , S(2).subs{:});
        end
    else
        error('dprior::subsref: unknown method (%s).', S(1).subs)
    end
  otherwise
    error('dprior::subsref: %s indexing not implemented.', S(1).type)
end
