function eqname = get_equation_name_by_number(eqnumber, M_)
% Returns the name of an equation given its number

% Copyright © 2022 Dynare Team
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

idx = find((cell2mat(M_.equations_tags(:,1)) == eqnumber) & ...
           strcmp(M_.equations_tags(:,2), 'name'));

if isempty(idx)
    eqname = '';
else
    eqname = M_.equations_tags{idx, 3};
end

end
