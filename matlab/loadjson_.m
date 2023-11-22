function o = loadjson_(jsonfilename)

% Reads a JSON file using jsondecode builtin.
%
% Returns the output in a format consistent with the output of the JSONlab
% toolbox, that we used to rely on (jsondecode builtin was introduced in MATLAB
% R2016b and in Octave 7). Now we no longer depend on JSONlab, but for
% historical reasons we still use its format.
%
% INPUTS
% - jsonfilename   [char]      1×n char array, name of the JSON file.
%
% OUTPUTS
% - o              [struct]    content of the JSON file.

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

json = fileread(jsonfilename);

o = jsondecode(json);

clear json

fnames = fieldnames(o);

for i=1:length(fnames)
    tmp = o.(fnames{i});
    if isstruct(tmp) && length(tmp)>1
        TMP = cell(length(tmp), 1);
        for j=1:length(tmp)
            TMP{j} = tmp(j);
        end
        o.(fnames{i}) = TMP;
    end
end
