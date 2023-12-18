function smoother(parameters, varlist)

% Computes smoothed variables.
%
% INPUTS
% - parameters    [char]    If row char array, possible values are 'posterior mode', 'posterior mean',
%                           'posterior median', 'mle mode', 'prior mode', 'prior mean' or 'calibration'.
% - varlist       [cell]    list of endogenous variables.
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

global M_ estim_params_ oo_ options_ bayestopt_

if nargin<2
    varlist = {};
end

parameters = strrep(parameters, ' ', '_');

[oo_, M_, options_, bayestopt_] = evaluate_smoother(parameters, varlist, M_, oo_, options_, bayestopt_, estim_params_);
