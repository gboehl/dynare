function [indices] = select_from_table(table,key,value)
% Copyright © 2010-2017 Dynare Team
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
candidates = table(strmatch(key,table(:,2),'exact'),:);
if nargin == 2
    indices = cell2mat( candidates(:,1) );
    return
end
indices = candidates(strmatch(value, candidates(:,3), 'exact'),1);
indices = cell2mat(indices);
