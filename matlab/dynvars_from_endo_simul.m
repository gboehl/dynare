function y2 = dynvars_from_endo_simul(y, it_, M_)
% Given the matrix y of paths for all endogenous (same format as
% oo_.endo_simul), and an iteration number (first simulation period corresponds
% to it_=M_.maximum_lag+1), return a vector of endogenous values in the format
% expected by the dynamic.m file (i.e. whose indices are described by
% M_.lead_lag_incidence)

% Copyright (C) 2020 Dynare Team
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

y2 = y(:, it_+(-M_.maximum_endo_lag:M_.maximum_endo_lead));
y2 = y2(find(M_.lead_lag_incidence'));
