function s=get_moments_size(options)
% function s=get_moments_size(options)
% Computes size of moments for a given model and options
%
% INPUTS
%    options: structure of Dynare options
%
% OUTPUTS
%    s: size of moments for a given model and options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2008-2017 Dynare Team
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

n = size(options.varlist,1);

if n == 0
    n = M_.endo_nbr;
end

n2 = n*n;

s = n; % mean
s = s + n;  % std errors
s = s + n2; % variance
s = s + n2; % correlations
s = s + options.ar*n2; % auto-correlations