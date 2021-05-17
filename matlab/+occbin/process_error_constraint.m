function error_parsed = process_error_constraint(constraint)
% function constraint1 = process_error_constraint(constraint)
% Constructs a string with the constraint error, i.e. by how much it is violated
% INPUTS
% - constraint      [char]     constraint to be parsed
%
% OUTPUTS
% - error_parsed    [char]     parsed constraint 

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

error_parsed = constraint;
error_parsed=regexprep(error_parsed,'=','');
num = length(regexp(error_parsed,'>'))+length(regexp(error_parsed,'<'));
error_parsed=regexprep(error_parsed,'>','-(');
error_parsed=regexprep(error_parsed,'<','-(');
for k=1:num
    error_parsed=[error_parsed ')'];
end

