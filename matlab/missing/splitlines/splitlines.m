function s = splitlines(string)

% SPLITLINES splits strings at newline characters into a string array.
%
% INPUT
% - string   [string]  String to be splitted with newline characters.
%
% OUTPUT
% - s   [string]
%
% Copyright © 2019 Dynare Team
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

assert(ischar(string) && ndims(string)==2 && size(string,1)<=1, 'The first argument has to be a row char array!');

s = strsplit(string, '\n');
s = reshape(s,size(s, 2), size(s,1));
end