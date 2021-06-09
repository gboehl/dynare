function o = convertjsondecode(o)

% Converts the output of jsondecode to be consistent with the output of loadjson.
%
% INPUTS
% - o     [struct]      Output of jsondecode.
%
% OUTPUTS
% - o     [struct]      Converted output of jsondecode.
%
% REMARKS
% The fields returned by the loadjson are systematically cell
% array, while jsondecode returns structure arrays if
% possible. This routine reorganize the data consistently with
% loadjson.

% Copyright (C) 2020 Dynare Team
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