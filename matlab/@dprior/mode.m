function m = mode(o, resetmoments)

% Return the prior mode.
%
% INPUTS
% - o               [dprior]
% - resetmoments    [logical]     Force the computation of the prior mode
%
% OUTPUTS
% - m               [double]      n×1 vector, prior mode

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

if nargin<2
    resetmoments = false;
end

if any(isnan(o.p5))
    resetmoments = true;
end

if resetmoments
    o.p5 = NaN(size(o.p5));
    o.moments('mode');
end

m = o.p5;
