% Stub implementation of datetime.
% It just returns a string, which is the equivalent of what one would
% get under MATLAB ≥ R2014b and a European locale with:
%  sprintf('%s', datetime)
%
% Hence, it only works as a substitute for the real datetime in the context of
% formatted output.

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

function d = datetime()
    d = datestr(now, 'dd-mmm-yyyy HH:MM:SS');
end
