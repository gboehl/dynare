function steady(printflag)

% Computes and prints the steady state.
%
% INPUTS
% - printflag    [logical]    scalar, print steady state if true (default value is true).
%
% OUTPUTS
% None

% Copyright © 2021 Dynare Team
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

global options_

if ~nargin || isempty(printflag)
    printflag = true;
end

noprint = options_.noprint;

options_.noprint = ~printflag;

steady();

options_.noprint = noprint;