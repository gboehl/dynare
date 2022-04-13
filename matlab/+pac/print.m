function print(pacexpectationmodelname, eqname, withcalibration)

% Prints the exansion of the PAC_EXPECTATION term in files.
%
% INPUTS
% - pacexpectationmodelname       [string]     Name of the expectation model.
% - eqname                         [string]    Name of the equation.
% - withcalibration                [logical]   Prints calibration if true.
%
% OUTPUTS
% None
%
% REMARKS
% The routine creates two text files
%
% - {pacexpectationmodelname}-parameters.inc     which contains the declaration of the parameters specific to the expectation model kind term.
% - {pacexpectationmodelname}-expression.inc     which contains the expanded version of the expectation model kind term.
%
% These routines are saved under the {modfilename}/model/pacexpectationmodel subfolder, and can be
% used after in another mod file (ie included with the macro directive @#include).

% Copyright © 2018-2019 Dynare Team
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

if nargin<3
    % Print calibration by default.
    withcalibration = true;
end

print_expectations(eqname, pacexpectationmodelname, 'pac', withcalibration);