function [irf_matching_file_name, irf_matching_file_path] = check_irf_matching_file(irf_matching_file)
% [irf_matching_file_name, irf_matching_file_path] = check_irf_matching_file(irf_matching_file)
% -------------------------------------------------------------------------
% Check if the provided irf_matching_file is a valid MATLAB function with
% .m extension and return name, path and extension of the file.
% -------------------------------------------------------------------------
% INPUTS
% - irf_matching_file:  [string] user provided name (with possible path and extension)
%                                of the MATLAB function that transforms model IRFs
% -------------------------------------------------------------------------
% OUTPUTS
% - irf_matching_file_name: [string] name of the MATLAB function (without extension)
% - irf_matching_file_path: [string] path to irf_matching_file_name
% -------------------------------------------------------------------------
% This function is called by
% - mom.run
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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


if isempty(irf_matching_file)
    % no irf_matching_file provided, so no transformations will be done
    irf_matching_file_name = '';
    irf_matching_file_path = '';
else
    [irf_matching_file_path, irf_matching_file_name, irf_matching_file_ext] = fileparts(irf_matching_file);
    % make sure file is a MATLAB function with .m extension
    if ~strcmp(irf_matching_file_ext,'.m')
        if strcmp(irf_matching_file_ext,'')
            irf_matching_file_ext = '.m';
        else
            error('method_of_moments: ''irf_matching_file'' needs to point towards a MATLAB function with extension ''.m''!');
        end
    end
    if isempty(irf_matching_file_path)
        irf_matching_file_path = '.';
    end
    if exist([irf_matching_file_path filesep irf_matching_file_name irf_matching_file_ext],'file') ~= 2
        error('method_of_moments: Could not find a ''irf_matching_file'' called ''%s''!',[irf_matching_file_path filesep irf_matching_file_name irf_matching_file_ext]);
    end
end