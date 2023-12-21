function y0 = get_mean(varargin)
% function y0 = get_mean(varargin)
% returns the mean of a variable identified by its name
%
% INPUTS:
%   vargargin           inputs containing
%                       - vname1, vname2, ... :  list of variable names
%                       - order: if integer 1 or 2, optionally last input can trigger the order
%                           at which steady state is computed
%
% OUTPUTS
%   y0:      mean values
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019-2023 Dynare Team
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

global M_ options_ oo_

y0=get_mean_no_globals(M_, oo_, options_, varargin{:});