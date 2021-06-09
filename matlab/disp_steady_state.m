function disp_steady_state(M,oo,options)
% function disp_steady_state(M,oo,options)
% computes and prints the steady state calculations
%
% INPUTS
%   M      structure of parameters
%   oo     structure of results
%   options structure of options
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2020 Dynare Team
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

skipline()
if options.loglinear
    disp('STEADY-STATE RESULTS FOR THE UNLOGGED VARIABLES:')
else
    disp('STEADY-STATE RESULTS:')
end
skipline()
endo_names = char(M.endo_names);
steady_state = oo.steady_state;

for i = 1:M.orig_endo_nbr
    fprintf('%s \t\t %g\n', endo_names(i,:), steady_state(i));
end
