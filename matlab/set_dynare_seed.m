function set_dynare_seed(varargin)
% set_dynare_seed(varargin)
% Set seeds depending on Matlab (Octave) version. This routine is
% a wrapper for set_dynare_seed_local_options

% Copyright © 2010-2023 Dynare Team
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

if nargin<1
    error('set_dynare_seed:: I need at least one input argument!')
end

options_=set_dynare_seed_local_options(options_,varargin{:});