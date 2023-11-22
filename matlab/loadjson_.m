function o = loadjson_(jsonfilename)

% Reads a json file using jsonlab toolbox or jsondecode builtin if available.
%
% INPUTS
% - jsonfilename   [char]      1×n char array, name of the JSON file.
%
% OUTPUTS
% - o              [struct]    content of the JSON file.
%
% REMARKS
% jsondecode builtin was introduced in MATLAB R2016b and in Octave 7.

% Copyright © 2020-2023 Dynare Team
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

if isoctave && octave_ver_less_than('7')
    o = loadjson(jsonfilename);
    return
end

json = fileread(jsonfilename);

o = jsondecode(json); clear('json');
o = convertjsondecode(o);
