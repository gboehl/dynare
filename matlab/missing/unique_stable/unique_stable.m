% Crude implementation of unique(…, 'stable'), which is missing in Octave < 6

% Copyright (C) 2021 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

function [r, idx] = unique_stable(x)

idx = [];

if iscell(x)
    r = {};
    for i = 1:numel(x)
        if ~ismember(x{i}, r)
            r{end+1} = x{i};
            idx(end+1) = i;
        end
    end
else
    r = [];
    for i = 1:numel(x)
        if ~ismember(x(i), r)
            r(end+1) = x(i);
            idx(end+1) = i;
        end
    end
end

if size(x, 2) == 1
    r = r';
end
