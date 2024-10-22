function mega = size_of_the_reduced_form_model(dr)
% Computes the size of dr.

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

names = fieldnames(dr);
number_of_scalars = 0;
for field=1:size(names,1)
    number_of_scalars = number_of_scalars + prod(size(getfield(dr,names{field})));
end
mega = 8 * number_of_scalars / 1048576 ;