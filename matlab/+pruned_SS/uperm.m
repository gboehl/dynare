function p = uperm(a)
% Return all unique permutations of possibly-repeating array elements
% =========================================================================
% Copyright © 2014 Bruno Luong <brunoluong@yahoo.com>
% Copyright © 2020 Dynare Team
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
% =========================================================================
% Original author: Bruno Luong <brunoluong@yahoo.com>, April 20, 2014
% https://groups.google.com/d/msg/comp.soft-sys.matlab/yQKVPTYrv6Q/gw1MzNd9sYkJ
% https://stackoverflow.com/a/42810388

[u, ~, J] = unique(a);
p = u(up(J, length(a)));

function p = up(J, n)
ktab = histc(J,1:max(J));
l = n;
p = zeros(1, n);
s = 1;
for i=1:length(ktab)
    k = ktab(i);
    c = nchoosek(1:l, k);
    m = size(c,1);
    [t, ~] = find(~p.');
    t = reshape(t, [], s);
    c = t(c,:)';
    s = s*m;
    r = repmat((1:s)',[1 k]);
    q = accumarray([r(:) c(:)], i, [s n]);
    p = repmat(p, [m 1]) + q;
    l = l - k;
end
end

end % uperm