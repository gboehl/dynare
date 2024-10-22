function hms = dynsec2hms(secs)
% DYNSEC2HMS Converts a number of seconds into a hours-minutes-seconds string

% Copyright © 2008-2009 Dynare Team
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

secs = round(secs);
s = rem(secs, 60);
m = rem(floor(secs / 60), 60);
h = floor(secs / 3600);
hms = sprintf('%dh%02dm%02ds', h, m, s);