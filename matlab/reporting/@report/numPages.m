function n = numPages(o)
%function n = numPages(o)
% return the number of pages currently in the report
%
% INPUTS
%   o     [report]  report object
%
% OUTPUTS
%   n     [integer] number of pages in the report object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2013-2015 Dynare Team
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

n = length(o.pages);
end