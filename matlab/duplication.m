function [Dp,DpMPinv] = duplication(p)
% [Dp,DpMPinv] = duplication(p)
% -------------------------------------------------------------------------
% Duplication Matrix as defined by Magnus and Neudecker (2002), p.49
% =========================================================================
% INPUTS
%   p:  [integer] length of vector
% -------------------------------------------------------------------------
% OUTPUTS
%   Dp:      Duplication matrix
%   DpMPinv: Moore-Penroze inverse of Dp
% -------------------------------------------------------------------------
% This function is called by
%   * identification.get_jacobians.m (previously getJJ.m)
% =========================================================================
% Copyright © 1997 Tom Minka <minka@microsoft.com>
% Copyright © 2019 Dynare Team
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
% Original author: Thomas P Minka (tpminka@media.mit.edu), April 22, 2013

a = tril(ones(p));
i = find(a);
a(i) = 1:length(i);
a = a + transpose(tril(a,-1));

j = a(:);

m = p*(p+1)/2;
Dp = spalloc(p*p,m,p^2);
for r = 1:size(Dp,1)
    Dp(r, j(r)) = 1;
end

if nargout > 1
    DpMPinv = (transpose(Dp)*Dp)\transpose(Dp);
end
