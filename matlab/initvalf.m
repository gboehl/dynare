function [series, p] = initvalf(DynareModel, options)

% Handles options for histvalf() and initvalf()
%
% INPUTS
% - caller           [char]      row array, name of calling function
% - DynareModel      [struct]    model description, a.k.a M_
% - options          [struct]    options specific to initivalf
%
% OUTPUTS
% - series           [dseries]   selected data from a file or a dseries
% - p                [integer]   number of periods (excluding the initial and terminal conditions)

% Copyright © 2003-2022 Dynare Team
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

[series, p] = histvalf_initvalf('INITVALF', DynareModel, options);