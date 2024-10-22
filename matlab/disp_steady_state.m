function disp_steady_state(M_,oo_,options_)
% function disp_steady_state(M_,oo_,options_)
% computes and prints the steady state calculations
%
% INPUTS
%   M_       structure of parameters
%   oo_      structure of results
%   options_ structure of options
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2001-2023 Dynare Team
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
if options_.loglinear
    disp('STEADY-STATE RESULTS FOR THE UNLOGGED VARIABLES:')
else
    disp('STEADY-STATE RESULTS:')
end
skipline()
endo_names = char(M_.endo_names);
steady_state = oo_.steady_state;

for i = 1:M_.orig_endo_nbr
    fprintf('%s \t\t %g\n', endo_names(i,:), steady_state(i));
end
