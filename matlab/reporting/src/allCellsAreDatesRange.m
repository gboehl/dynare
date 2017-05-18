function tf = allCellsAreDatesRange(dcell)
%function tf = allCellsAreDatesRange(dcell)
% Determines if all the elements of dcell are a range of dates
%
% INPUTS
%   dcell     cell of dates
%
% OUTPUTS
%   tf        true if every entry of dcell is a range of dates
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2014-2015 Dynare Team
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

assert(iscell(dcell));
tf = true;
for i=1:length(dcell)
    if ~(isdates(dcell{i}) && dcell{i}.ndat >= 2)
        tf = false;
        return;
    end
end
end