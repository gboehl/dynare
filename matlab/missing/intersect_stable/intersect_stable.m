% Crude implementation of intersect(…, 'stable'), which is missing in Octave

% Copyright (C) 2019 Dynare Team
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

function [c, ia, ib] = intersect_stable(a, b)

if nargin ~= 2
    error('intersect_stable: needs exactly 2 input arguments');
end

if isnumeric (a) && isnumeric (b)
    c = [];
elseif iscell (b)
    c = {};
else
    c = '';
end
ia = [];
ib = [];

if isempty (a) || isempty (b)
    return
else
    isrowvec = isrow (a) && isrow (b);

    for i = 1:numel(a)
        if iscellstr(c)
            idx = strcmp(a(i), b);
        else
            idx = a(i) == b;
        end
        if any(idx) && ~ismember(a(i), c)
            c = [c(:); a(i)];
            if nargout > 1
                ia = [ia, i];
                ib = [ib, find(idx)];
            end
        end
    end

    %% Adjust output orientation for MATLAB compatibility
    if isrowvec
        c = c.';
    end
end
end

%!test
%! a = [3 4 1 5];
%! b = [2 4 9 1 6];
%! [c,ia,ib]=intersect_stable(a,b);
%! assert(c, [4 1])
%! assert(ia, [2 3])
%! assert(ib, [2 4])
%! assert(a(ia), c)
%! assert(b(ib), c)

%!test
%! a = [3 4 1 5]';
%! b = [2 4 9 1 6]';
%! [c,ia,ib]=intersect_stable(a,b);
%! assert(c, [4 1]')
%! assert(ia, [2 3])
%! assert(ib, [2 4])
%! assert(a(ia), c)
%! assert(b(ib), c)

%!test
%! a = { 'defun', 'mapcar', 'let', 'eval-when'};
%! b = { 'setf', 'let', 'list', 'cdr', 'defun'};
%! [c,ia,ib]=intersect_stable(a,b);
%! assert(c, { 'defun', 'let' })
%! assert(ia, [1 3])
%! assert(ib, [5 2])
%! assert(a(ia), c)
%! assert(b(ib), c)
